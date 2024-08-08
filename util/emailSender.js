// util/emailSender.js

const nodemailer = require('nodemailer');
const ejs = require('ejs');
const path = require('path');
const archiver = require('archiver');
const fs = require('fs');
require('dotenv').config(); // 确保加载 .env 文件

const sendEmail = async (to) => {
    let subject="xx网站的结果"
    const templateData = {
        title: "您好！，您的数据已经运行完成",
        message: "您的数据结果在压缩包内"
    };
    try {
        // 指定要发送的文件
        const files = ['img_1.png'];

        // 创建压缩包并将文件添加到压缩包中
        const zipFileName = `${new Date().toISOString().split('T')[0]}.zip`;
        const zipFilePath = path.join(__dirname, '../', zipFileName);

        const archive = archiver('zip', {
            zlib: { level: 9 } // 设置压缩等级
        });

        const output = fs.createWriteStream(zipFilePath);

        archive.pipe(output);

        files.forEach(file => {
            const filePath = path.join(__dirname, '../public/images', file);
            if (fs.existsSync(filePath)) {
                archive.file(filePath, { name: file });
            } else {
                console.warn(`文件 ${filePath} 不存在，跳过。`);
            }
        });

        await archive.finalize();

        // 配置邮件发送器
        let transporter = nodemailer.createTransport({
            host: 'smtp.vip.163.com',
            port: 465,
            secure: true, // 使用SSL
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
        try {
            let info = await transporter.sendMail(mailOptions);
            console.log('邮件发送成功: ' + info.response);
            // 删除压缩包文件
            fs.unlinkSync(zipFilePath);
        } catch (error) {
            console.error('发送邮件时出错: ' + error);
            // 删除压缩包文件
            fs.unlinkSync(zipFilePath);
            throw error; // 重新抛出错误以便外部捕获
        }
    } catch (error) {
        console.error('处理邮件时出错: ' + error);
        throw error; // 重新抛出错误以便外部捕获
    }
};

module.exports = sendEmail;
