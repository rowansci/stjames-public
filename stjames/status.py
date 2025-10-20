from enum import Enum


class Status(int, Enum):
    # what a job gets when it's created if user is below max_concurrency
    QUEUED = 0

    # still running...
    RUNNING = 1

    # what a job gets when it has finished
    COMPLETED_OK = 2

    # something went wrong
    FAILED = 3

    # job was stopped externally, maybe timed out or something
    STOPPED = 4

    # Status if a user has exceeded their allowed max_concurrency
    AWAITING_QUEUE = 5

    # Job is not yet submitted
    DRAFT = 6
