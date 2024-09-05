// util/pythonScriptRunner.js
const { exec } = require('child_process');
const sendEmail = require('./emailSender');
const db = require('../config/db'); // 确保已连接到数据库
const fs = require('fs-extra');
const path = require('path');

/**
 * 运行Python脚本的函数
 * @param {number} type - 表示要运行的Python脚本类型
 * @param {string} folderPath - Python脚本使用的文件夹路径参数
 * @param {string} storagePath - 存储路径参数
 * @param {string} email - 用于发送邮件的电子邮件地址
 * @param {number} scriptRunId - 数据库中对应脚本运行的ID
 * @returns {Promise<string>} - 返回一个Promise，解析为Python脚本输出的生成目录路径
 */
const runPythonScript = (fun,type, folderPath, storagePath, email, scriptRunId) => {
    console.log("正在执行py脚本，请等待");
    return new Promise((resolve, reject) => {
        // 脚本映射
        const scriptMap = {
            1: `./py/${fun}/singleCell.py`,
            2: `./py/${fun}/singleCellSpatial.py`,
            3: `./py/${fun}/BTSpatial.py`,
            4: `./py/${fun}/Xenium.py`,
            5: `./py/${fun}/h5ad.py`
        };

        // 根据type获取对应的Python脚本路径
        const scriptPath = scriptMap[type];
        if (!scriptPath) {
            reject(`未知的type: ${type}，无法找到对应的Python脚本`);
            return;
        }

        // 指定Python解释器路径
        const pythonInterpreter = process.env.PYTHON_INTERPRETER;
        // 创建运行Python脚本的命令
        const pythonCommand = `${pythonInterpreter} ${scriptPath} ${folderPath} ${storagePath}`;
        exec(pythonCommand, async (error, stdout, stderr) => {
            if (stderr.includes("UserWarning: No data for colormapping provided") ||
                stderr.includes("TypeError: close() argument must be a Figure")) {
            } else if (error){
                console.error(`执行Python脚本出错: ${error.message}`);
                await updateScriptRunStatus(scriptRunId, 'failed'); // 更新状态为失败
                reject(`执行Python脚本出错: ${stderr}`);
                return;
            }
            console.log(`Python脚本over`)
            const generatedDirPath = stdout.trim();
            console.log(`生成的目录路径: ${generatedDirPath}`)

            try {
                // 发送邮件
                console.log(`发送邮件...`)
                await sendEmail(email, generatedDirPath);
                console.log('邮件发送成功');

                // 更新数据库中的状态为 'completed'
                await updateScriptRunStatus(scriptRunId, 'completed');
                console.log('数据库状态更新为 completed');
                const folderToDelete1 = path.dirname(storagePath);
                const folderToDelete2 = folderToDelete1.replace('storage', 'uploads');

                await fs.rm(folderToDelete1, { recursive: true, force: true });
                await fs.rm(folderToDelete2, { recursive: true, force: true });
                resolve(generatedDirPath);
            } catch (emailError) {
                console.error('发送邮件时出错:', emailError);
                await updateScriptRunStatus(scriptRunId, 'failed'); // 更新状态为失败
                reject(`发送邮件失败: ${emailError.message}`);
            }
        });
    });
};
const updateScriptRunStatus = async (id, status) => {
    await db.query(`
        UPDATE script_runs SET status = ? WHERE id = ?
    `, [status, id]);
};

module.exports = { runPythonScript };
