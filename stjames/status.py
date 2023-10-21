from enum import Enum


class Status(int, Enum):
    # what a job gets when it's created
    QUEUED = 0

    # still running...
    RUNNING = 1

    # what a job gets when it has finished
    COMPLETED_OK = 2

    # something went wrong
    FAILED = 3

    # job was stopped externally, maybe timed out or something
    STOPPED = 4
