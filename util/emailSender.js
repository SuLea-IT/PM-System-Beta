// util/emailSender.js

const nodemailer = require('nodemailer');
const ejs = require('ejs');
const path = require('path');

const sendEmail = async (to, subject, templateData) => {
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
                path: path.join(__dirname, '../public/images/img_1.png')
            }
        ]
    };

    // 发送邮件
    try {
        let info = await transporter.sendMail(mailOptions);
        console.log('邮件发送成功: ' + info.response);
    } catch (error) {
        console.error('发送邮件时出错: ' + error);
    }
};

module.exports = sendEmail;
