import sys

from config import rfam_local as config


def main(db_name):
    db = getattr(config, db_name.upper())
    options = "--host {host} --port {port} --user {user} -p{password} --database {db_name}".format(
        host=db['host'],
        port=db['port'],
        user=db['user'],
        password=db['pwd'],
        db_name=db['db'],
    )
    sys.stdout.write(options)


if __name__ == '__main__':
    main(sys.argv[1])
