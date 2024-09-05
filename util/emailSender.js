// util/emailSender.js

const nodemailer = require('nodemailer');
const ejs = require('ejs');
const path = require('path');
const archiver = require('archiver');
const fs = require('fs');
require('dotenv').config();

const sendEmail = async (to, dirPath) => {
    let subject = "xx网站的结果";
    const templateData = {
        title: "您好！，您的数据已经运行完成",
        message: "您的数据结果在压缩包内"
    };

    try {
        const zipFileName = `${new Date().toISOString().split('T')[0]}.zip`;
        const zipFilePath = path.join(__dirname, '../', zipFileName);

        const archive = archiver('zip', {
            zlib: { level: 9 }
        });

        const output = fs.createWriteStream(zipFilePath);

        archive.pipe(output);
        console.log(dirPath)
        // 使用完整路径压缩指定的文件夹内的所有文件和子文件夹
        archive.directory(path.resolve(dirPath), false);

        await new Promise((resolve, reject) => {
            output.on('close', resolve);
            archive.on('error', reject);
            archive.finalize();
        });

        // 配置邮件发送器
        let transporter = nodemailer.createTransport({
            host: 'smtp.vip.163.com',
            port: 465,
            secure: true,
            auth: {
                user: process.env.EMAIL_USER,
                pass: process.env.EMAIL_PASS
            }
        });

        // 渲染模板
        const templatePath = path.join(__dirname, 'templates', 'emailTemplate.ejs');
        const htmlContent = await ejs.renderFile(templatePath, templateData);

        // 邮件选项
        let mailOptions = {
            from: process.env.EMAIL_USER,
            to: to,
            subject: subject,
            html: htmlContent,
            attachments: [
                {
                    filename: zipFileName,
                    path: zipFilePath
                }
            ]
        };

        // 发送邮件
        let info = await transporter.sendMail(mailOptions);
        console.log('邮件发送成功: ' + info.response);
        fs.unlinkSync(zipFilePath); // 删除压缩包文件
    } catch (error) {
        console.error('处理邮件时出错: ' + error);
        throw error;
    }
};

module.exports = sendEmail;
