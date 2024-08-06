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

// For tracking upload progress
const uploadProgress = {};
const fileExistFlag = {}; // To flag if the file already exists

// Calculate MD5 for a file chunk
const calculateMD5 = (data) => {
    return crypto.createHash('md5').update(data).digest('hex');
};

// Handle file upload and merging
async function handleFileUpload(req) {
    const {index, totalChunks, fileName} = req.body;
    const total = parseInt(totalChunks, 10);
    if (isNaN(total) || total < 1) {
        throw new Error('Invalid total chunks number');
    }

    const projectId = req.body.projectId;
    const userId = req.body.userId;
    const uploadId = `${projectId}_${userId}_${fileName}`;

    const projectDir = `uploads/${projectId}`;
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
        return {code: 400, msg: 'File already exists', data: null};
    }

    if (!req.file) {
        throw new Error('No file provided');
    }

    try {
        const fileChunkSize = req.file.size;
        uploadProgress[uploadId].fileSize += fileChunkSize;

        const chunkData = await fsExtra.readFile(req.file.path);
        const chunkMD5 = calculateMD5(chunkData);

        if (parseInt(index, 10) === 0) {
            const [existingFiles] = await db.query('SELECT id FROM files WHERE first_chunk_md5 = ?', [chunkMD5]);
            if (existingFiles.length > 0) {
                fileExistFlag[uploadId] = true;
                await fsExtra.remove(chunkDir);
                await fsExtra.remove(req.file.path);
                return {code: 400, msg: 'File already exists', data: null};
            }
            uploadProgress[uploadId].firstChunkMD5 = chunkMD5;
        }

        if (await fsExtra.pathExists(chunkPath)) {
            return {code: 400, msg: 'Chunk already exists', data: null};
        }

        await fsExtra.move(req.file.path, chunkPath);
        const indexInt = parseInt(index, 10);
        if (isNaN(indexInt) || indexInt < 0 || indexInt >= uploadProgress[uploadId].chunks.length) {
            throw new Error('Invalid index value');
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
                return {code: 400, msg: 'File already exists', data: null};
            }

            await fsExtra.move(tempOutputPath, finalOutputPath);

            await db.query('INSERT INTO files (project_id, user_id, filename, filepath, file_size, first_chunk_md5) VALUES (?, ?, ?, ?, ?, ?)', [
                projectId,
                userId,
                fileName,
                finalOutputPath,
                uploadProgress[uploadId].fileSize,
                uploadProgress[uploadId].firstChunkMD5
            ]);

            await fsExtra.remove(chunkDir);

            delete uploadProgress[uploadId];
            delete fileExistFlag[uploadId];

            return {code: 200, msg: 'File successfully merged and chunk directory deleted', data: {md5: finalMD5}};
        } else {
            return {code: 200, msg: 'Chunk uploaded successfully, waiting for other chunks', data: null};
        }
    } catch (error) {
        if (error.message.includes('dest already exists')) {
            await fsExtra.remove(chunkDir);
            await fsExtra.remove(req.file.path);
            return {code: 400, msg: 'File already exists', data: null};
        } else {
            await fsExtra.remove(chunkDir);
            throw new Error(`Server error: unable to move file chunk - ${error.message}`);
        }
    }
}

module.exports = {
    handleFileUpload
};
