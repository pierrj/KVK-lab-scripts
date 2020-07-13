import sys
import csv
import sqlite3

coverage_file = str(sys.argv[1])

scaffold_number = int(sys.argv[2])

with open(coverage_file) as coverage:
    coverage_reader = csv.reader(coverage, delimiter = '\t')
    for row in coverage_reader:
        for i in range(scaffold_number):
            conn = sqlite3.connect(r"scaffold"+str(i+1)+"_sql.db")
            c = conn.cursor()
            c.execute('''Drop TABLE if exists server''')
            c.execute('''Create TABLE if not exists server(base, count)''')
            to_add = []
            while int(row[0]) == i+1:
                to_add.append((int(float(row[1])), int(float(row[2]))))
                try:
                    row = next(coverage_reader)
                except StopIteration:
                    break
            c.executemany("INSERT INTO server(base, count) VALUES(?,?)", to_add)
            c.execute("CREATE INDEX base_index ON server(base)")
            conn.commit()
            conn.close()