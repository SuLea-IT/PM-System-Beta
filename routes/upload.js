//routes/upload.js
const express = require('express');
const multer = require('multer');
const path = require('path');
const {handleFileUpload} = require('../util/fileUploadHandler'); // 引入FileUpload模块

const router = express.Router();
const multerUpload = multer({dest: 'uploads/'});
const uploadGene = {
    /*
    * 1.代表单细胞数据类型
    * 2.单细胞级别空间类型
    * 3.百迈克空间转录组数据
    * 4.Xenium数据
    * 5.h5ad数据类型
    * */
    1: {
        allowedExtensions: ['.tsv.gz', '.mtx.gz','.txt','.text'],
        requiredFileNames: ['barcodes', 'features', 'matrix', "*"],
        uploadFileCount: 4
    },
    2: {
        allowedExtensions: ['.tsv.gz', '.mtx.gz','.txt','.text'],
        requiredFileNames: ['barcodes', 'features', 'matrix', 'barcodes_pos', "*"],
        uploadFileCount: 5
    },
    3: {
        allowedExtensions: ['.tsv.gz', '.mtx.gz','.txt','.text'],
        requiredFileNames: ['barcodes', 'features', 'matrix', '*'],
        uploadFileCount: 5
    },
    4: {
        allowedExtensions: ['.csv.gz','.h5','.txt','.text'],
        requiredFileNames: ['*', '*'],
        uploadFileCount: 3
    },
    5: {
        allowedExtensions: ['.h5ad','.txt','.text'],
        requiredFileNames: ['*', '*'],
        uploadFileCount: 2
    }
};
const uploadGenes = {
    /*
    * 1.代表单细胞数据类型
    * 2.单细胞级别空间类型
    * 3.百迈克空间转录组数据
    * 4.Xenium数据
    * 5.h5ad数据类型
    * */
    1: {
        allowedExtensions: ['.tsv.gz', '.mtx.gz','.txt','.text'],
        requiredFileNames: ['barcodes', 'features', 'matrix', "*"],
        uploadFileCount: 4
    },
    2: {
        allowedExtensions: ['.tsv.gz', '.mtx.gz','.txt','.text'],
        requiredFileNames: ['barcodes', 'features', 'matrix', 'barcodes_pos', "*"],
        uploadFileCount: 5
    },
    3: {
        allowedExtensions: ['.tsv.gz', '.mtx.gz','.txt','.text'],
        requiredFileNames: ['barcodes', 'features', 'matrix', '*'],
        uploadFileCount: 5
    },
    4: {
        allowedExtensions: ['.csv.gz','.h5','.txt','.text'],
        requiredFileNames: ['*', '*'],
        uploadFileCount: 3
    },
    5: {
        allowedExtensions: ['.h5ad','.txt','.text'],
        requiredFileNames: ['*', '*'],
        uploadFileCount: 2
    }
};
const uploadCluster = {
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

function validateFileTypeAndName(fun,fileName, type) {
    let restrictions;
    if(fun==1){
        restrictions = uploadCluster[type];
    }else if(fun==2){
        restrictions = uploadGene[type];
    }
    else if(fun==3){
        restrictions = uploadGenes[type];
    }

    // 获取完整的扩展名
    const fileExtension = fileName.slice(fileName.indexOf('.'));  // 获取第一个点之后的所有字符
    const fileNameWithoutExtension = fileName.slice(0, fileName.indexOf('.')).toLowerCase();

    if (!restrictions.allowedExtensions.includes(fileExtension)) {
        throw new Error(`文件类型不支持。只允许上传以下类型: ${restrictions.allowedExtensions.join(', ')}`);
    }

    let matchFound = false;
    for (const required of restrictions.requiredFileNames) {
        if (required === '*' || fileNameWithoutExtension.includes(required.toLowerCase())) {
            matchFound = true;
            break;
        }
    }

    if (!matchFound) {
        throw new Error(`文件名必须包含以下之一: ${restrictions.requiredFileNames.join(', ')}`);
    }

    return true;
}

router.post('/upload', multerUpload.array('file', 10), async (req, res) => {
    try {
        const {type, fileName,fun} = req.body;
        const uploadResults = [];
        let Fdata = "上传中";
        let FCode = 200;
        let restrictions
        if(fun==1){
            restrictions = uploadCluster[type];
        }else if(fun==2){
            restrictions = uploadGene[type];
        }else if(fun==3){
            restrictions = uploadGenes[type];
        }


        // 校验传递的 fileName 字段
        validateFileTypeAndName(fun,fileName, type);

        // 获取完整的文件扩展名
        const fileExtension = fileName.slice(fileName.indexOf('.'));  // 获取第一个点之后的所有字符

        for (const file of req.files) {
            req.file = file;
            const result = await handleFileUpload(req, fileExtension); // 传递fileExtension参数
            // console.log(result)

            if (result.msg === "文件合并成功并且分片文件夹已删除") {
                Fdata = "所有文件上传成功";
                uploadResults.push(result.data);
            } else if (result.msg === '文件已经存在') {
                Fdata = "文件已经存在";
                FCode = 400;
                break;
            }
        }

        res.status(FCode).json({
            code: FCode,
            msg: Fdata,
            data: uploadResults
        });
    } catch (error) {
        console.error('File upload error:', error);
        res.status(400).json({
            code: 400,
            msg: error.message,
            data: null
        });
    }
});


module.exports = router;
