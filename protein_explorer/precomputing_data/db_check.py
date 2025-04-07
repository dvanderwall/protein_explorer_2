import sqlite3

# Connect to the database
conn = sqlite3.connect("C:\\Users\\mz30\\AppData\\Local\\Temp\\blosum_analysis.db")
cursor = conn.cursor()

# Check if the table exists
cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='blosum_scores'")
tables = cursor.fetchall()
print(f"Tables found: {tables}")

# Check if there's data in the table
if tables:
    cursor.execute("SELECT COUNT(*) FROM blosum_scores")
    count = cursor.fetchone()[0]
    print(f"Total rows in blosum_scores: {count}")
    
    # If there's data, get a sample
    if count > 0:
        cursor.execute("SELECT score FROM blosum_scores LIMIT 10")
        sample = cursor.fetchall()
        print(f"Sample data: {sample}")
else:
    print("No blosum_scores table found")

conn.close()