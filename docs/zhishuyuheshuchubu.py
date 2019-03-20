import math

def get_zhishu(num):
    a = True
    n = int(math.sqrt(num))+1
    for i in range(2, n):
        if num % n == 0:
            print("%d is not a zhishu" % num)
            a = False
            break
    if a is True:
        print("%d is a zhishu" % num)

get_zhishu(1356)
