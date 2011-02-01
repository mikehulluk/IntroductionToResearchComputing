from sqlalchemy import *
from sqlalchemy.orm import *

# Create Database:
db = create_engine('sqlite:///tutorial.db')
metadata = MetaData(db)

users = Table('users', metadata,
    Column('user_id', Integer, primary_key=True),
    Column('name', String(40)),
    Column('age', Integer),
    Column('password', String),
)

users.create()

# Add some data:
i = users.insert()
i.execute(name='Mary', age=30, password='secret')
i.execute({'name': 'John', 'age': 42},
          {'name': 'Susan', 'age': 57},
          {'name': 'Carl', 'age': 33})

# Find some data:
s = users.select()
rs = s.execute()

for row in rs:
    print row.name, 'is', row.age, 'years old'
