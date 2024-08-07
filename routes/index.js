// routes/index.js

let express = require('express');
let router = express.Router();

const usersRouter = require('./users');
const uploadRouter = require('./upload');
const emailRouter = require('./email');
/* GET home page. */
router.get('/', function(req, res, next) {
  res.render('index', { title: 'Express' });
});


// 配置子路由
router.use('/users', usersRouter);
router.use('/upload', uploadRouter);
router.use('/email', emailRouter);
module.exports = router;

