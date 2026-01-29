from .base import Base, LowercaseStrEnum


class MessageType(LowercaseStrEnum):
    """Message severity level."""

    ERROR = "error"
    WARNING = "warning"
    INFO = "info"


class Message(Base):
    """User-facing message with severity, title, and body."""

    type: MessageType
    title: str
    body: str
