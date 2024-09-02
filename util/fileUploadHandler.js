//util\fileUploadHandler.js
const fsExtra = require('fs-extra');
const fs = require('fs-extra');
const path = require('path');
const db = require('../config/db');
const crypto = require('crypto');
const { exec } = require('child_process');
const sendEmail = require('./emailSender');
require('dotenv').config();
const {
    checkChunkExists,
    mergeChunks,
    createMD5Incremental,
    updateMD5Incremental,
    finalizeMD5Incremental
} = require('./fileHelpers');
const scriptMap = {
    1: './py/singleCell.py',
    2: './py/singleCellSpatial.py',
    3: './py/BTSpatial.py',
    4: './py/Xenium.py',
    5: './py/h5ad.py'
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
        allowedExtensions: ['.csv.gz','.h5'],
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
const checkAllFilesUploaded = (type,totalFiles, currentFileIndex) => {
    console.log("检测开始")
    // 检查所有文件是否都上传完成
    const allUploaded = Object.values(uploadProgress).every(progress => progress.chunks.every(status => status === true));
    console.log("检测完成")
    if (allUploaded) {
        if(currentFileIndex==uploadRestrictions[type].uploadFileCount){
            console.log("检测通过")

            return true;
        }
    }else {
        console.log("检测不通过")
    }
    return false; // 如果没有全部上传完成或条件不满足，返回false
};
const createStorageStructure = async (timestamp, ip, type1,type) => {
    console.log("开始", process.cwd())
    // 创建根目录下的 storage 文件夹
    const storageDir = path.join(process.cwd(), 'storage');
    await fs.ensureDir(storageDir);

    // 在 storage 文件夹下创建以当前日期为名称的子文件夹
    const dateDir = path.join(storageDir, timestamp);
    await fs.ensureDir(dateDir);

    // 在日期子文件夹下创建以 ip_小时.分钟 为名称的进一步子文件夹
    const finalDir = path.join(dateDir, `${ip}_${type}`);
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
const runPythonScript = (type, folderPath, storagePath) => {
    console.log("正在执行py脚本，请等待");
    return new Promise((resolve, reject) => {
        // 根据type获取对应的Python脚本路径
        const scriptPath = scriptMap[type];
        if (!scriptPath) {
            reject(`未知的type: ${type}，无法找到对应的Python脚本`);
            return;
        }

        // 指定Python解释器路径
        const pythonInterpreter = process.env.PYTHON_INTERPRETER;
        // 传递 storagePath 参数给 Python 脚本
        const pythonCommand = `${pythonInterpreter} ${scriptPath} ${folderPath} ${storagePath}`;

        exec(pythonCommand, (error, stdout, stderr) => {
            if (error) {
                console.error(`执行Python脚本出错: ${error.message}`);
                reject(`执行Python脚本出错: ${stderr}`);
                return;
            }
            // 捕获Python脚本的输出，并将其解析为文件路径
            const generatedDirPath = stdout.trim();
            resolve(generatedDirPath);
        });
    });
};


const handleFileUpload = async (req, fileExtension) => {
    const { index, totalChunks, fileName, type, number, currentFileIndex, totalFiles, email } = req.body;
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
    const uploadId = `${timestamp}_${ip}_${type}_${sanitizedFileName}`;

    const projectDir = path.join('uploads', timestamp);
    await fsExtra.ensureDir(projectDir);

    const chunkDir = path.join(projectDir, uploadId);
    await fsExtra.ensureDir(chunkDir);

    const chunkPath = path.join(chunkDir, indexInt.toString());

    if (!uploadProgress[uploadId]) {
        uploadProgress[uploadId] = {
            chunks: new Array(total).fill(false),
            fileName,
            fileSize: 0,
            scriptExecuted: false,
        };
        uploadProgress[uploadId].md5Incremental = createMD5Incremental();
        fileExistFlag[uploadId] = false;
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

        const allUploaded = uploadProgress[uploadId].chunks.every(status => status === true);
        if (allUploaded) {
            const files = uploadProgress[uploadId].chunks.map((_, i) => path.join(chunkDir, i.toString()));
            const tempOutputPath = path.join(projectDir, `temp_${uploadId}`);
            await mergeChunks(files, tempOutputPath);

            const finalMD5 = finalizeMD5Incremental(uploadProgress[uploadId].md5Incremental);

            const finalFolder = path.join(projectDir, `${ip}_${type}`);
            await fsExtra.ensureDir(finalFolder);

            const finalOutputPath = path.join(finalFolder, `${fileName}`);

            if (await fsExtra.pathExists(finalOutputPath)) {
                return { code: 400, msg: '目标文件已经存在', data: null };
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
                const storagePath = await createStorageStructure(timestamp, ip, type1, type);
                console.log(storagePath);

                setTimeout(async () => {
                    try {
                        const generatedFilePaths = await runPythonScript(type, finalFolder, storagePath);
                        console.log(`Python脚本生成的文件路径: ${generatedFilePaths}`);

                        await sendEmail(email, storagePath);

                        const folderToDelete1 = path.dirname(storagePath);
                        const folderToDelete2 = folderToDelete1.replace('storage', 'uploads');

                        await fs.rm(folderToDelete1, { recursive: true, force: true });
                        await fs.rm(folderToDelete2, { recursive: true, force: true });

                        // 清理缓存
                        delete uploadProgress[uploadId];
                        delete fileExistFlag[uploadId];

                    } catch (error) {
                        console.error('处理脚本或发送邮件时出错:', error);
                        // 在发送邮件或脚本执行出错时也清理缓存
                        delete uploadProgress[uploadId];
                        delete fileExistFlag[uploadId];
                    }
                }, 0);
            }

            return { code: 200, msg: '文件合并成功并且分片文件夹已删除', data: { md5: finalMD5, first_chunk_md5: firstChunkMD5 } };
        } else {
            return { code: 200, msg: '分片上传成功，等待其他分片', data: null };
        }
    } catch (error) {
        if (error.message.includes('dest already exists')) {
            console.error(`文件上传错误: ${error.message}`, error);
            await fsExtra.remove(chunkDir);
            await fsExtra.remove(req.file.path);
            return { code: 400, msg: '文件已经存在', data: null };
        } else {
            console.error(`移动文件分片错误: ${error.message}`, error);
            await fsExtra.remove(chunkDir);
            // 清理缓存
            delete uploadProgress[uploadId];
            delete fileExistFlag[uploadId];
            throw new Error(`服务器错误：无法移动文件分片 - ${error.message}`);
        }
    }
};




module.exports = {
    handleFileUpload
};
