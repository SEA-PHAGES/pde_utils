"""
Functions to multithread process a list of inputs.
"""

from queue import LifoQueue
import threading


def create_threads(target, work_stack, stack_lock, num_threads):
    class MixedThread(threading.Thread):
        def __init__(self, threadid):
            threading.Thread.__init__(self)

            self.name = threadid

            self.target = target
            self.stack = work_stack
            self.lock = stack_lock

        def run(self):
            while not self.stack.empty():
                self.lock.acquire()

                work_item = self.stack.get()
                self.lock.release()

                self.target(*work_item)

    threads = []
    for x in range(num_threads):
        threads.append(MixedThread(x+1))

    return threads


def multithread(target, work_items, threads=1):
    lock = threading.Lock()
    stack = LifoQueue()

    for item in work_items:
        stack.put(item)

    threads = create_threads(target, stack, lock, threads)

    for thread in threads:
        thread.start()

    for thread in threads:
        thread.join()
