# **Rfam Production Code**


### **About**
--------
This repository provides a set of tools related to the Rfam production pipeline.
The collection includes script for data export, database processing, visualization
and validation.

### **Scripts**
--------
* **Export:** Data Export tools
* **Processing:** Database processing tools
* **View:** Rfam family view process related scripts
* **Validation:** Data validation scripts
* **Utils:** Support modules/utilities

### Installation

```
# install Python dependencies
pip install -r requirements.txt

# set up Python path
export PYTHONPATH=/path/to/the/project

# create and edit configuration file (excluded from version control)
cp config/rfam_local_template.py config/rfam_local.py
```

### **Requirements**
---------------
* Python 2.6 or later
* Json 1.1.1
* Python Mysql Connector 2.1.3
* Python MySQL 0.7.1

### **To Do**
--------
* Clan Competition validation script
* View process validation script
* `utils/RfamDB.py` to be converted to a class

### **Contact Us**
-------------
Feel free to contact us via email rfam-help@ebi.ac.uk or submit an [issue](https://github.com/Rfam/rfam-production/issues).
