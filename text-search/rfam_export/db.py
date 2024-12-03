# -*- coding: utf-8 -*-

# Copyright [2009-2024] EMBL-European Bioinformatics Institute
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


@frozen
class RfamDb:
    mysql_host: str
    mysql_user: str
    mysql_password: str
    mysql_port: int
    mysql_database: str

    @contextmanager
    def connection(self):
        connection = pymysql.connect(
            host=self.mysql_host,
            user=self.mysql_user,
            password=self.mysql_password,
            port=self.mysql_port,
            database=self.mysql_database,
            cursorclass=pymysql.cursors.DictCursor,
        )
        yield connection
        connection.close()

    @contextmanager
    def cursor(self):
        with self.connection() as conn:
            with conn.cursor() as cursor:
                yield cursor
