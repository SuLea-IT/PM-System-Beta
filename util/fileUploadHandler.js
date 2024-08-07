const fsExtra = require('fs-extra');
const path = require('path');
const db = require('../config/db');
const crypto = require('crypto');
const {
    checkChunkExists,
    mergeChunks,
    createMD5Incremental,
    updateMD5Incremental,
    finalizeMD5Incremental
} = require('./fileHelpers');

const uploadProgress = {};
const fileExistFlag = {};

const calculateMD5 = (data) => {
    return crypto.createHash('md5').update(data).digest('hex');
};

const formatDate = (date) => {
    const year = date.getFullYear();
    const month = date.getMonth() + 1;
    const day = date.getDate();

    return `${year}-${month}-${day}`;
};

const sanitizeIP = (ip) => {
    if (ip.includes(':')) {
        // 检查是否是IPv6的本地回环地址
        if (ip === '::1') {
            return '127.0.0.1';
        }
        // 处理其他IPv6地址
        try {
            const ip6To4 = require('ip6-to4');
            const ipv4 = ip6To4(ip);
            if (ipv4) {
                return ipv4;
            } else {
                // 如果转换失败，保留原始IPv6地址
                return ip;
            }
        } catch (error) {
            // 如果ip6-to4转换失败，保留原始IPv6地址
            return ip;
        }
    }
    return ip; // 返回原始的IPv4地址
};


async function handleFileUpload(req) {
    const {index, totalChunks, fileName} = req.body;
    const total = parseInt(totalChunks, 10);
    if (isNaN(total) || total < 1) {
        throw new Error('无效的分片总数');
    }

    const timestamp = formatDate(new Date());

    let ip = req.headers['x-forwarded-for'] || req.connection.remoteAddress;
    ip = sanitizeIP(ip);

    const uploadId = timestamp + '_' + ip + '_' + fileName;

    const projectDir = `uploads/${timestamp}`;
    await fsExtra.ensureDir(projectDir);

    const chunkDir = path.join(projectDir, uploadId);
    const chunkPath = path.join(chunkDir, index.toString());
    await fsExtra.ensureDir(chunkDir);

    if (!uploadProgress[uploadId]) {
        uploadProgress[uploadId] = {chunks: new Array(total).fill(false), fileName, fileSize: 0};
        uploadProgress[uploadId].md5Incremental = createMD5Incremental();
        fileExistFlag[uploadId] = false;
    }

    if (fileExistFlag[uploadId]) {
        await fsExtra.remove(chunkDir);
        return {code: 400, msg: '文件已经存在', data: null};
    }

    if (!req.file) {
        throw new Error('文件未提供');
    }

    try {
        const fileChunkSize = req.file.size;
        uploadProgress[uploadId].fileSize += fileChunkSize;

        const chunkData = await fsExtra.readFile(req.file.path);
        const chunkMD5 = calculateMD5(chunkData);

        if (parseInt(index, 10) === 0) {
            console.log('第一个分片的MD5:', chunkMD5);
            const [existingFiles] = await db.query('SELECT id FROM files WHERE first_chunk_md5 = ?', [chunkMD5]);
            if (existingFiles) {
                console.log('文件已经存在，阻止上传');
                fileExistFlag[uploadId] = true;
                await fsExtra.remove(chunkDir);
                await fsExtra.remove(req.file.path);
                return {code: 400, msg: '文件已经存在', data: null};
            }
            uploadProgress[uploadId].firstChunkMD5 = chunkMD5;
        }

        if (await fsExtra.pathExists(chunkPath)) {
            return {code: 400, msg: '分片已经存在', data: null};
        }

        await fsExtra.move(req.file.path, chunkPath);
        const indexInt = parseInt(index, 10);
        if (isNaN(indexInt) || indexInt < 0 || indexInt >= uploadProgress[uploadId].chunks.length) {
            throw new Error('无效的索引值');
        }

        updateMD5Incremental(uploadProgress[uploadId].md5Incremental, chunkData);

        uploadProgress[uploadId].chunks[indexInt] = true;

        const allUploaded = uploadProgress[uploadId].chunks.every(status => status === true);
        if (allUploaded) {
            const files = uploadProgress[uploadId].chunks.map((_, i) => path.join(chunkDir, i.toString()));
            const tempOutputPath = path.join(projectDir, `temp_${uploadId}`);
            await mergeChunks(files, tempOutputPath);

            const finalMD5 = finalizeMD5Incremental(uploadProgress[uploadId].md5Incremental);
            const finalOutputPath = path.join(projectDir, `${finalMD5}${path.extname(fileName)}`);

            if (await fsExtra.pathExists(finalOutputPath)) {
                await fsExtra.remove(chunkDir);
                await fsExtra.remove(finalOutputPath);
                return {code: 400, msg: '目标文件已经存在', data: null};
            }

            await fsExtra.move(tempOutputPath, finalOutputPath);

            await db.query('INSERT INTO files (timestamp, ip, filename, filepath, file_size, first_chunk_md5) VALUES (?, ?, ?, ?, ?, ?)', [
                timestamp,
                ip,
                fileName,
                finalOutputPath,
                uploadProgress[uploadId].fileSize,
                uploadProgress[uploadId].firstChunkMD5
            ]);

            await fsExtra.remove(chunkDir);

            delete uploadProgress[uploadId];
            delete fileExistFlag[uploadId];

            return {code: 200, msg: '文件合并成功并且分片文件夹已删除', data: {md5: finalMD5}};
        } else {
            return {code: 200, msg: '分片上传成功，等待其他分片', data: null};
        }
    } catch (error) {
        if (error.message.includes('dest already exists')) {
            console.error(`文件上传错误: ${error.message}`, error);
            await fsExtra.remove(chunkDir);
            await fsExtra.remove(req.file.path);
            return {code: 400, msg: '文件已经存在', data: null};
        } else {
            console.error(`移动文件分片错误: ${error.message}`, error);
            await fsExtra.remove(chunkDir);
            throw new Error(`服务器错误：无法移动文件分片 - ${error.message}`);
        }
    }
}

module.exports = {
    handleFileUpload
};
