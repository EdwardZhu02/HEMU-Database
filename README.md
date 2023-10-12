# The HEMU Andropogoneae Database and Analysis Platform

## 1 Project Introduction

HEMU is an easy-to-handle, graphical based database and analysis platform for Androposoneae grasses which enables user to discover genomic data from multiple representative Andropogoneae species within a click.

The name "HEMU" draws its inspiration from the Chinese classic, the Book of Songs, symbolizing robust crops and bountiful harvests. HEMU encompasses massive multi-omics data and is equipped with six sophisticated analysis toolkits, which enables easy utilization of the bioinformatics data and customized performance of comparative analysis from novel perspectives among Andropogoneae species.

The project is currently deployed and available for open research use [HERE](https://shijunpenglab.com/HEMUdb/).  
To cite this project: [HEMU: an integrated Andropogoneae comparative genomics database and analysis platform | bioRxiv](https://www.biorxiv.org/content/10.1101/2023.05.19.541421v1)

## 2 Deploying the Main HEMU Framework

Currently, the project is deployed on a server running Ubuntu 20.04 LTS.

**Step1**: create environment and install prerequisites

```bash
conda env create -f HEMUdb_new.yaml
```

**Step2**: configure MySQL server, creating a database for data storage

```mysql
CREATE DATABASE hemu_database DEFAULT CHARSET utf8mb4 COLLATE utf8mb4_general_ci;
```

**Step3**: configure Redis server, which can be run with docker

```bash
docker run --name redis-hemu -p 6379:6379 redis
```

**Step4**: Start the server instance using uwsgi, then start the celery task backend

```bash
# Start server using uwsgi configurations
uwsgi --ini /path/to/config/hemu_uwsgi.ini

# Start celery worker for handling async tasks
cd /path/to/projectdir/ && nohup celery -A HEMU_Database_Main worker -l info > /path/to/logdir/hemu_celery.log 
```
