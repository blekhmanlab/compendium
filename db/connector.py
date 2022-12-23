import sqlite3
import time

class Connection(object):
    """Data type holding the data required to maintain a database
    connection and perform queries.

    """
    def __init__(self):
        """Stores db connection info in memory and initiates a
        connection to the specified db."""
        try:
            self.db = sqlite3.connect("compendium.db")
        except Exception as e:
            print(f'FATAL: {e}')
            exit(1)
        print('Connected!')
        self.setup_tables()

    def write(self, query, params=None):
        cursor = self.db.cursor()
        if params is not None:
            cursor.execute(query, params)
        else:
            cursor.execute(query)
        self.db.commit()
        cursor.close()

    def read(self, query, params=None):
        """Helper function that converts results returned stored in a
        sqlite3 cursor into a less temperamental list format.

        Arguments:
            - query: The SQL query to be executed.
            - params: Any parameters to be substituted into the query.
                sqlite3 handles this better than Python does.
        Returns:
            - A list of tuples, one for each row of results.

        """

        results = []
        try:
            cursor = self.db.cursor()
            if params is not None:
                cursor.execute(query, params)
            else:
                cursor.execute(query)

            for result in cursor:
                results.append(result)
            
            cursor.close()

            return results

        except Exception as e:
            print(f'ERROR with db query execution: {e}')

    def setup_tables(self):
        self.write("""
            CREATE TABLE IF NOT EXISTS samples(
                srs TEXT PRIMARY KEY,
                project TEXT,
                taxon TEXT,
                srr TEXT,
                library_strategy TEXT,
                library_source TEXT,
                pubdate TEXT,
                total_bases INTEGER,
                exported INTEGER DEFAULT 0
            )
        """)

        self.write("""
            CREATE TABLE IF NOT EXISTS tags (
                tagid INTEGER PRIMARY KEY,
                srs TEXT,
                tag TEXT,
                value TEXT
            )
        """)

        self.write("""
            CREATE TABLE IF NOT EXISTS tags (
                tagid INTEGER PRIMARY KEY,
                srs TEXT,
                tag TEXT,
                value TEXT
            )
        """)

    def __del__(self):
        """Closes the database connection when the Connection object
        is destroyed."""

        if self.db is not None:
            self.db.close()