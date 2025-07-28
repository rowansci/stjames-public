from pydantic import BaseModel


class SCFSettings(BaseModel):
    max_iters: int = 250

    soscf: bool = False
