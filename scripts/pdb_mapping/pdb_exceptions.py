class PDBException(Exception):
    pass


class QCChecksError(PDBException):
    def __init__(self):
        message = "Error when running QC checks on PDB information in the database"
        super().__init__(message)


class CheckRowsError(QCChecksError):
    def __init__(self, num_rows_pdb_temp, num_rows_pdb):
        message = "There is too large a difference in the number of rows in the table." \
                  " num_rows_pdb_temp: {0} num_rows_pdb: {1}".format(num_rows_pdb_temp, num_rows_pdb)
        super().__init__(message)


class NoFileGivenError(PDBException):
    def __init__(self):
        message = "Please provide a text file with the information to import to the database."
        super().__init__(message)
