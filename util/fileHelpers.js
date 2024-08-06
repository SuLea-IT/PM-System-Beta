const crypto = require('crypto');
const fsExtra = require('fs-extra');
const path = require('path');

// Initialize an incremental MD5 object
function createMD5Incremental() {
    return crypto.createHash('md5');
}

// Update incremental MD5 object
function updateMD5Incremental(md5Incremental, chunk) {
    md5Incremental.update(chunk, 'utf8');
}

// Get the final MD5 value
function finalizeMD5Incremental(md5Incremental) {
    return md5Incremental.digest('hex');
}

// Check if a chunk exists
async function checkChunkExists(chunkDir, chunkFilename) {
    const fullPath = path.join(chunkDir, chunkFilename);
    try {
        await fsExtra.access(fullPath);
        return true;
    } catch (error) {
        return false;
    }
}

// Merge file chunks and calculate the final MD5
async function mergeChunks(files, dest) {
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
