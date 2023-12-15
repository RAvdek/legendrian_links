import time
from threading import Thread


class TimeoutSignaller(Thread):
    def __init__(self, limit_secs, handler):
        Thread.__init__(self)
        self.limit_secs = limit_secs
        self.running = True
        self.handler = handler
        assert callable(handler), "Timeout Handler needs to be a method"

    def run(self):
        timeout_limit = time.time() + self.limit_secs
        while self.running:
            if time.time() >= timeout_limit:
                self.handler()
                self.stop_run()
                break

    def stop_run(self):
        self.running = False


class ProcessContextManager:
    def __init__(self, process, seconds=0):
        self.seconds = seconds
        self.process = process
        self.signal = TimeoutSignaller(self.seconds, self.signal_handler)

    def __enter__(self):
        self.signal.start()
        return self.process

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.signal.stop_run()

    def signal_handler(self):
        # Make process terminate however you like
        # using self.process reference
        raise TimeoutError("Process took too long to execute")


import multiprocessing as mp

class RunWithTimer(object):

    def __init__(self, fun, args, t_secs):
        self.results = []
        self.fun = fun
        self.args = args
        self.t_secs = t_secs

    def f_mod(self, args):
        print(args)
        self.results.append(self.fun(args))

    def run(self):
        t_proc = mp.Process(target=time.sleep, args=(self.t_secs,))
        f_proc = mp.Process(target=self.f_mod, args=(self.args,))
        t_proc.start()
        f_proc.start()
        while f_proc.is_alive():
            if not t_proc.is_alive():
                print("TIMEOUT")
                f_proc.terminate()
                t_proc.terminate()
            print(self.results)
        return self.results

def returnx(x):
    time.sleep(.1)
    return x

if __name__ == "__main__":
    print(RunWithTimer(returnx, (1,), 10).run())

