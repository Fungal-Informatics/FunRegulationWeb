import schedule
import time

def job():
    print("job: I'm working...")


def setup():
    print("cron_setup: setting up cron")
    schedule.every(5).seconds.do(job)

    while True:
        schedule.run_pending()
        time.sleep(1)
