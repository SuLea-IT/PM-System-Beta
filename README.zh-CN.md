# 项目管理系统测试版
简体中文 | [English](https://github.com/SuLea-IT/PM-System-Beta/blob/main/README.md)

> 这个项目是一个基于 `Node` / `Python` 的后端应用，主要用于部分数据分析。

## 目录结构

```bash
PM-System-Beta/
├── PM-System-Beta.iml
├── README.md
├── README.zh-CN.md
├── app.js
├── bin/
│   └── www
├── package-lock.json
├── package.json
├── public/
│   └── stylesheets/
├── routes/
│   ├── index.js
│   └── users.js
├── sql/
│   └── PM-System.sql
└── views/
    ├── error.pug
    ├── index.pug
    └── layout.pug

```

1. 克隆项目代码

   ```cmd
    git clone https://github.com/SuLea-IT/PM-System-Beta.git
   ```

2. 进入项目目录

   ```cmd
   cd PM-System-Beta
   ```

3. 安装依赖

   ```
   npm install
   ```

4. 运行项目

   ```cmd
   node bin/www
   ```
## 拉黑机制
### 1. IP访问计数

设某IP地址的访问次数为 \( N_{ip} \)，则：

\[N_{ip} = \text{访问次数}\]

### 2. IP拉黑逻辑

设某IP地址的拉黑状态为 \( B_{ip} \)，如果该IP地址在一天内的访问次数 \( N_{ip} \) 超过300次，则：

\[B_{ip} = \begin{cases}1 & \text{若 } N_{ip} > 300 \\0 & \text{其他}\end{cases}\]

拉黑时间 \( T_{ip} \) 为1小时（3600秒）：

\[T_{ip} = 3600 \quad \text{秒}\]

### 3. 国家IP拉黑逻辑

设某国家的被拉黑的IP地址列表为 \( L_c \)，并设该国家在一个月内被拉黑的IP地址数量为 \( |L_c| \)。如果 \( |L_c| \) 超过20，则该国家标记为受到限制：

\[|L_c| > 20\]

对于受到限制的国家，其所有IP地址的访问限制变为250次，设 \( N_{ip,\text{limit}} \) 为访问限制次数：

\[N_{ip,\text{limit}} = \begin{cases}250 & \text{若 } |L_c| > 20 \\300 & \text{其他}\end{cases}\]

### 4. 国家拉黑升级逻辑

设同一国家的IP地址被拉黑的次数为 \( N_{c, \text{blacklist}} \)，若超过一定次数，则增加拉黑时间。

拉黑时间 $$\( T_{ip, \text{upgrade}} \) $$递增公式为：
$$

\[T_{ip, \text{upgrade}} = T_{ip} + 600 \times (N_{c, \text{blacklist}} - \text{threshold})\]

$$
其中，`threshold` 是定义的基础拉黑次数阈值。




## 功能流程图：
- 大文件分片上传(增量计算)：

![img_2.png](public/images/img_2.png)
