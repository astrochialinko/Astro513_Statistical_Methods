
import time

class TimerError(Exception):
    """An exception used to report errors in use of Timer class"""

class Timer:

    def __init__(self, text="", name=None):
        self.start_time = None
        self.name = name
        self.elapsed_time = 0

    def clear(self):
       self.start_time = None
       self.elapsed_time = 0

    def start(self):
        if self.start_time is not None:
            raise TimerError(f"Timer[{self.name}].start: Timer is running.")
        elif self.start_time is None:
            self.start_time = time.perf_counter()

    def cstart(self):
        self.clear()
        self.start()

    def stop(self):
        if self.start_time is None:
            raise TimeoutError(f"Timer[{self.name}].stop: Timer is not running.")

        # Stop the timer then report the elapsed time
        elif self.start_time is not None:
            elapsed_time = time.perf_counter() - self.start_time
            self.start_time = None
            self.elapsed_time += elapsed_time

    def report(self, text):
        print(f"{text} {self.elapsed_time}")

    def elapsed(self):
        return self.elapsed_time

if __name__ == "__main__":

    t1 = Timer(name="t1")


    t1.cstart()
    time.sleep(1)
    t1.stop()
    time.sleep(1)
    t1.start()
    time.sleep(1)
    t1.stop()

    t1.report("total time:")
