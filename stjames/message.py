from .base import Base, LowercaseStrEnum


class MessageType(LowercaseStrEnum):
    ERROR = "error"
    WARNING = "warning"
    INFO = "info"


class Message(Base):
    type: MessageType
    title: str
    body: str
