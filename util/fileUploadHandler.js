//util\fileUploadHandler.js
const fsExtra = require('fs-extra');
const fs = require('fs-extra');
const path = require('path');
const db = require('../config/db');
const crypto = require('crypto');
const {exec} = require('child_process');
require('dotenv').config();
const { runPythonScript } = require('./pythonScriptRunner');
const {
    checkChunkExists,
    mergeChunks,
    createMD5Incremental,
    updateMD5Incremental,
    finalizeMD5Incremental
} = require('./fileHelpers');
const saveScriptRun = async (type,fun, folderPath, storagePath, email) => {
    try {
        const result = await db.query(`
            INSERT INTO script_runs (type,fun, folderPath, storagePath, email)
            VALUES (?,?, ?, ?, ?)
        `, [type,fun, folderPath, storagePath, email]);
        return result.insertId; // 返回新插入的记录ID
    } catch (error) {
        console.error('插入脚本运行参数时出错:', error.message);
        throw error; // 抛出错误，供上层处理
    }
};

const uploadRestrictions = {
    /*
    * 1.代表单细胞数据类型
    * 2.单细胞级别空间类型
    * 3.百迈克空间转录组数据
    * 4.Xenium数据
    * 5.h5ad数据类型
    * */
    1: {
        allowedExtensions: ['.tsv.gz', '.mtx.gz'],
        requiredFileNames: ['barcodes', 'features', 'matrix'],
        uploadFileCount: 3
    },
    2: {
        allowedExtensions: ['.tsv.gz', '.mtx.gz'],
        requiredFileNames: ['barcodes', 'features', 'matrix', 'barcodes_pos'],
        uploadFileCount: 4
    },
    3: {
        allowedExtensions: ['.tsv.gz', '.mtx.gz'],
        requiredFileNames: ['barcodes', 'features', 'matrix', '*'],
        uploadFileCount: 4
    },
    4: {
        allowedExtensions: ['.csv.gz', '.h5'],
        requiredFileNames: ['*', '*'],
        uploadFileCount: 2
    },
    5: {
        allowedExtensions: ['.h5ad'],
        requiredFileNames: ['*', '*'],
        uploadFileCount: 1
    }
};
const uploadProgress = {};
const fileExistFlag = {};

const calculateMD5 = (data) => {
    return crypto.createHash('md5').update(data).digest('hex');
};
const checkAllFilesUploaded = (type, totalFiles, currentFileIndex) => {
    console.log("检测开始")
    // 检查所有文件是否都上传完成
    if (currentFileIndex == totalFiles) {
        console.log("检测通过")

        return true;
    } else {
        console.log("检测不通过")
    }
    return false; // 如果没有全部上传完成或条件不满足，返回false
};
const getScriptRunById = async (id) => {
    try {
        const rows = await db.query(`SELECT * FROM script_runs WHERE id = ?`, [id]);
        if (rows.length > 0) {
            return rows[0];
        } else {
            console.log(`未找到脚本运行ID: ${id}`);
            return null; // 或者抛出错误
        }
    } catch (error) {
        console.error('读取脚本运行记录时出错:', error.message);
        throw error; // 抛出错误，供上层处理
    }
};

const createStorageStructure = async (timestamp, ip, type1, type,fun) => {
    console.log("开始", process.cwd())
    // 创建根目录下的 storage 文件夹
    const storageDir = path.join(process.cwd(), 'storage');
    await fs.ensureDir(storageDir);

    // 在 storage 文件夹下创建以当前日期为名称的子文件夹
    const dateDir = path.join(storageDir, timestamp);
    await fs.ensureDir(dateDir);

    // 在日期子文件夹下创建以 ip_小时.分钟 为名称的进一步子文件夹
    const finalDir = path.join(dateDir, `${ip}_${fun}_${type}`);
    await fs.ensureDir(finalDir);

    const zfinalDir = path.join(finalDir, `${type}`);
    await fs.ensureDir(zfinalDir);

    // 获取相对于项目根目录的路径
    const relativePath = path.relative(process.cwd(), zfinalDir);
    console.log("最后", relativePath)
    return relativePath;
};
const formatDate = (date) => {
    const year = date.getFullYear();
    const month = date.getMonth() + 1;
    const day = date.getDate();

    return `${year}-${month}-${day}`;
};

const sanitizeIP = (ip) => {
    if (ip.includes(':')) {
        if (ip === '::1') {
            return '127.0.0.1';
        }

        // Handle IPv6-mapped IPv4 addresses (e.g., ::ffff:192.168.1.1)
        if (ip.startsWith('::ffff:')) {
            return ip.split(':').pop(); // Extract the IPv4 part
        }

        try {
            const ip6To4 = require('ip6-to4');
            const ipv4 = ip6To4(ip);
            return ipv4 || ip;
        } catch (error) {
            return ip;
        }
    }
    return ip;
};


const getCurrentHourAndMinute = (date) => {
    const hours = date.getHours().toString().padStart(2, '0');
    const minutes = date.getMinutes().toString().padStart(2, '0');
    return `${hours}.${minutes}`;
};
const sanitizeFileName = (name) => {
    return name.replace(/[^a-zA-Z0-9-_]/g, '_');
};


const handleFileUpload = async (req, fileExtension) => {
    const {index, totalChunks, fileName, type, number, currentFileIndex, totalFiles, email,fun} = req.body;
    const type1 = type;
    const total = parseInt(totalChunks, 10);
    const indexInt = parseInt(index, 10);

    if (isNaN(total) || total < 1) {
        throw new Error('无效的分片总数');
    }

    if (isNaN(indexInt) || indexInt < 0 || indexInt >= total) {
        throw new Error('无效的索引值');
    }

    const timestamp = formatDate(new Date());
    let ip = req.headers['x-forwarded-for'] || req.connection.remoteAddress;
    ip = sanitizeIP(ip);

    const sanitizedFileName = sanitizeFileName(fileName);
    const uploadId = `${timestamp}_${ip}_${fun}_${type}_${sanitizedFileName}`;

    const projectDir = path.join('uploads', timestamp);
    await fsExtra.ensureDir(projectDir);
    await cleanupPreviousUpload(uploadId, projectDir);
    const chunkDir = path.join(projectDir, uploadId);
    await fsExtra.ensureDir(chunkDir);

    const chunkPath = path.join(chunkDir, indexInt.toString());

    if (!uploadProgress[uploadId]) {
        uploadProgress[uploadId] = {
            chunks: new Array(total).fill(false),
            fileName,
            fileSize: 0,
            scriptExecuted: false,
            md5Incremental: createMD5Incremental() // 在初始化时创建一个新的MD5对象
        };
        fileExistFlag[uploadId] = false;
    } else if (indexInt === 0) {
        // 如果是新的一次上传，重新初始化md5对象
        uploadProgress[uploadId].md5Incremental = createMD5Incremental();
    }


    if (!req.file) {
        throw new Error('文件未提供');
    }

    try {
        const fileChunkSize = req.file.size;
        uploadProgress[uploadId].fileSize += fileChunkSize;

        const chunkData = await fsExtra.readFile(req.file.path);
        const chunkMD5 = calculateMD5(chunkData);

        if (indexInt === 0) {
            console.log('第一个分片的MD5:', chunkMD5);
            uploadProgress[uploadId].firstChunkMD5 = chunkMD5;
        }

        await fsExtra.move(req.file.path, chunkPath);

        updateMD5Incremental(uploadProgress[uploadId].md5Incremental, chunkData);
        uploadProgress[uploadId].chunks[indexInt] = true;
        console.log("创建合并目录了")
        const allUploaded = uploadProgress[uploadId].chunks.every(status => status === true);
        if (allUploaded) {
            console.log(chunkDir)
            console.log(projectDir)
            await fsExtra.ensureDir(chunkDir);
            const files = uploadProgress[uploadId].chunks.map((_, i) => path.join(chunkDir, i.toString()));
            const tempOutputPath = path.join(projectDir, `temp_${uploadId}`);
            await mergeChunks(files, tempOutputPath);
            console.log("合并文件结束")
            const finalMD5 = finalizeMD5Incremental(uploadProgress[uploadId].md5Incremental);

            const finalFolder = path.join(projectDir, `${ip}_${fun}_${type}`);
            await fsExtra.ensureDir(finalFolder);

            const finalOutputPath = path.join(finalFolder, `${fileName}`);

            if (await fsExtra.pathExists(finalOutputPath)) {
                return {code: 400, msg: '目标文件已经存在', data: null};
            }

            await fsExtra.move(tempOutputPath, finalOutputPath);

            await db.query('INSERT INTO files (timestamp, ip, filename, filepath, file_size, first_chunk_md5, request_type, file_count) VALUES (?, ?, ?, ?, ?, ?, ?, ?)', [
                timestamp,
                ip,
                fileName,
                finalOutputPath,
                uploadProgress[uploadId].fileSize,
                uploadProgress[uploadId].firstChunkMD5,
                type,
                number,
            ]);

            await fsExtra.remove(chunkDir);

            const firstChunkMD5 = uploadProgress[uploadId].firstChunkMD5;
            if (checkAllFilesUploaded(type, totalFiles, currentFileIndex)) {
                const storagePath = await createStorageStructure(timestamp, ip, type1, type,fun);
                console.log(`创建的存储路径: ${storagePath}`);

                const scriptRunId = await saveScriptRun(type,fun, finalFolder, storagePath, email);
                console.log(`保存的脚本运行ID: ${scriptRunId}`);

                // 在尝试读取之前，等待插入完全完成
                await new Promise(resolve => setTimeout(resolve, 100)); // 增加一个短暂的延时，等待数据库事务完成

                // 从数据库读取脚本运行参数
                const scriptRun = await getScriptRunById(scriptRunId);

                // 检查并运行脚本
                if (scriptRun) {
                    runPythonScript(scriptRun.fun,scriptRun.type, scriptRun.folderPath, scriptRun.storagePath,
                    scriptRun.email, scriptRunId)
                        .then((result) => console.log(`Python脚本执行成功: ${result}`))
                        .catch((error) => console.error(`Python脚本执行失败: ${error}`));
                } else {
                    console.error('脚本运行记录未找到，可能在保存时出现问题。');
                }
            }

            return {code: 200, msg: '文件合并成功并且分片文件夹已删除', data: {md5: finalMD5, first_chunk_md5: firstChunkMD5}};
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
            delete uploadProgress[uploadId];
            delete fileExistFlag[uploadId];
            throw new Error(`服务器错误：无法移动文件分片 - ${error.message}`);
        }
    }
};


module.exports = {
    handleFileUpload
};
