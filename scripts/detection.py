from typing import Callable
import re

def is_aba_free(x: str) -> bool:
    return re.search(r"(.+)(.+)\1", x) is None

def is_abacaba_free(x: str) -> bool:
    return re.search(r"(.+)(.+)\1(.+)\1\2\1", x) is None

def is_square_free(x: str) -> bool:
    return re.search(r"(.+)\1", x) is None

def find_abacaba(x: str):
    matches = re.finditer(r"(.+)(.+)\1(.+)\1\2\1", x)
    for match in matches:
        yield match

def get_search_protocol(protocol: str) -> Callable:
    return {
            "abacaba": is_abacaba_free,
            "aba": is_aba_free,
            "square": is_square_free
            }[protocol]
