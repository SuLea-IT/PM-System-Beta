CREATE DATABASE DataManagement;

USE DataManagement;

CREATE TABLE DataDetails (
                             id INT AUTO_INCREMENT PRIMARY KEY,
                             category VARCHAR(255) NOT NULL,
                             name VARCHAR(255) NOT NULL,
                             level VARCHAR(255),
                             storage_path VARCHAR(255),
                             type ENUM('0', '1', '2', '3', '4', '5') NOT NULL
);

CREATE TABLE UploadRecords (
                               id INT AUTO_INCREMENT PRIMARY KEY,
                               visitor_ip VARCHAR(45) NOT NULL,
                               data_name VARCHAR(255) NOT NULL,
                               data_size INT NOT NULL,
                               upload_time TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
