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

CREATE TABLE files (
                       id INT AUTO_INCREMENT PRIMARY KEY,
                       project_id VARCHAR(255) NOT NULL,
                       user_id VARCHAR(255) NOT NULL,
                       filename VARCHAR(255) NOT NULL,
                       filepath VARCHAR(255) NOT NULL,
                       file_size INT NOT NULL,
                       first_chunk_md5 VARCHAR(255) NOT NULL
);
