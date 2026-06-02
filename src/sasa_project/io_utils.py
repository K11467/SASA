import csv


def read_csv_dicts(path):
    """Read CSV files written on Linux/Windows with a small encoding fallback."""
    last_error = None
    for encoding in ("utf-8", "gbk"):
        try:
            with open(path, encoding=encoding, newline="") as f:
                return list(csv.DictReader(f))
        except UnicodeDecodeError as exc:
            last_error = exc
    raise last_error
