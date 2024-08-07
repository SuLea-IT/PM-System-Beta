// routes/emailRoute.js

const express = require('express');
const router = express.Router();
const sendEmail = require('../util/emailSender');

router.post('/send-email', async (req, res) => {
    const { to, subject, title, message } = req.body;

    const templateData = {
        title: title,
        message: message
    };

    try {
        await sendEmail(to, subject, templateData);
        res.status(200).send('邮件发送成功！');
    } catch (error) {
        res.status(500).send('发送邮件时出错。');
    }
});

module.exports = router;
