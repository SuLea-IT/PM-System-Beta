// routes/emailRoute.js

const express = require('express');
const router = express.Router();
const sendEmail = require('../util/emailSender');

router.post('/send-email', async (req, res) => {
    const { to} = req.body;

    try {
        await sendEmail(to);
        res.status(200).send('邮件发送成功！');
    } catch (error) {
        console.error('发送邮件时出错: ', error);
        res.status(500).send(`发送邮件时出错: ${error.message}`);
    }
});

module.exports = router;
