const express = require('express');
const multer = require('multer');
const { handleFileUpload } = require('../util/fileUploadHandler'); // 引入FileUpload模块

const router = express.Router();
const upload = multer({ dest: 'temp/' }); // Set up multer for file handling
const multerUpload = multer({dest: 'uploads/'});
// 上传文件到项目
router.post('/upload',multerUpload.array('file', 10), async (req, res) => {
    try {
        const uploadResults = [];
        let Fdata = "所有文件上传成功";
        let FCode = 200;

        for (const file of req.files) {
            req.file = file;
            const result = await handleFileUpload(req);

            if (result.msg === "分片上传成功，等待其他分片") {
                uploadResults.push(result.data);
            } else if (result.msg === '文件已经存在') {
                Fdata = "文件已经存在";
                FCode = 400;
                break;  // 终止后续文件上传
            }
        }

        res.status(FCode).json({
            code: FCode,
            msg: Fdata,
            data: uploadResults
        });
    } catch (error) {
        console.error('File upload error:', error);
        res.status(500).json({
            code: 500,
            msg: '文件上传失败',
            data: error.message
        });
    }
});



module.exports = router;
