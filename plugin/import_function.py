import time
import datetime
import inspect
from cachetools import cached, TTLCache  # 1 - let's import the "cached" decorator and the "TTLCache" object from cachetools
cache = TTLCache(maxsize=100, ttl=300)  # 2 - let's create the cache object.
import os.path
import csv


#@cached(cache)  # 3 - it's time to decorate the method to use our cache system!
def get_candy_price(candy_id):
    # [1] /Users/Ethan/PycharmProjects/Python/Practice/Cache_test.py
    # [2] /Users/Ethan/PycharmProjects/Python/Practice/
    # [3] Cache_test.csv
    # [4] /Users/Ethan/PycharmProjects/Python/Practice/data/
    # [5] /Users/Ethan/PycharmProjects/Python/Practice/data/tmp_Cache_test.csv
    location = inspect.stack()[1][1]  # [1]
    directory = os.getcwd() + '/'  # [2]
    filename = 'tmp_' + location.replace(directory, "").replace(".py", ".csv")  # [3]
    save_directory = directory + 'data/'  # [4]
    save_file = save_directory + filename  # [5]
    #print(save_file)
    result = 123
    inputdata = (1, 2, 3e20)
    func = 'multiplication_factor'
    with open(save_file, newline='', encoding='utf-8-sig') as file2:
        rows = csv.DictReader(file2)
        for index, row in enumerate(rows):
            print(index, row['function_name'], row['result'])
    with open(save_file, 'w') as file:
        SaveData = csv.writer(file, delimiter=',')
        SaveData.writerow(['function_name', 'result', 'input'])
        SaveData.writerow([func, result, inputdata])
        print('test')

        #print(inspect.getmodule(inspect.stack()[1][0]).__file__)
    # let's use a sleep to simulate the time your function spends trying to connect to
    # the web service, 5 seconds will be enough.
    time.sleep(5)
    # let's pretend that the price returned by the web service is $1 for candies with a
    # odd candy_id and $1,5 for candies with a even candy_id
    price = 1.5 if candy_id % 2 == 0 else 1

    return (datetime.datetime.now().strftime("%c"), price)