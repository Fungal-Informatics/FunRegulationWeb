import database
import cron

print("BEGIN")
# exit(1)


database.setup()
cron.setup()