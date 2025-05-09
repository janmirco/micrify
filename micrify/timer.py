import logging
import time


class Timer:
    def __init__(self, name):
        self.name = name

    def __enter__(self):
        self.start = time.perf_counter()
        return self

    def __exit__(self, *_):
        self.end = time.perf_counter()

    def print_time(self, milli: bool = True):
        if milli:
            self.time_factor = 1000.0
            self.units = "milliseconds"
        else:
            self.time_factor = 1.0
            self.units = "seconds"
        self.time = self.time_factor * (self.end - self.start)
        logging.info(f"{self.name} took {self.time:.6f} {self.units}.")
