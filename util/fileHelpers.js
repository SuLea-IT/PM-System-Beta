const crypto = require('crypto');
const fsExtra = require('fs-extra');
const path = require('path');
let isFinalized = false;

// 初始化一个增量计算 MD5 的对象
const createMD5Incremental = () => {
    isFinalized = false; // 每次重新创建哈希对象时重置标志
    return crypto.createHash('md5');
};


// 更新增量 MD5 对象
function updateMD5Incremental(md5Incremental, chunk) {
    if (isFinalized) {
        throw new Error('Cannot update MD5 after finalization');
    }
    md5Incremental.update(chunk, 'utf8');
}


// 获取增量 MD5 的最终值
function finalizeMD5Incremental(md5Incremental) {
    if (!isFinalized) {
        isFinalized = true;  // 设置为已完成
        return md5Incremental.digest('hex');
    } else {
        return md5Incremental.digest('hex');
    }
}

// 检查文件分片是否存在
async function checkChunkExists(chunkDir, chunkFilename) {
    const fullPath = path.join(chunkDir, chunkFilename);
    try {
        await fsExtra.access(fullPath);
        return true;
    } catch (error) {
        return false;
    }
}

// 合并文件分片并计算整体文件的 MD5
async function mergeChunks(files, dest) {
    console.log("开始合并文件分片")
    return new Promise((resolve, reject) => {
        const output = fsExtra.createWriteStream(dest);

        function appendFile(file, callback) {
            const input = fsExtra.createReadStream(file);
            input.pipe(output, {end: false});
            input.on('end', callback);
            input.on('error', callback);
        }

        (function next(i) {
            if (i < files.length) {
                appendFile(files[i], (err) => {
                    if (err) {
                        reject(err);
                    } else {
                        next(i + 1);
                    }
                });
            } else {
                output.end();
                resolve();
            }
        })(0);
    });
}

module.exports = {
    createMD5Incremental,
    updateMD5Incremental,
    finalizeMD5Incremental,
    checkChunkExists,
    mergeChunks
};
