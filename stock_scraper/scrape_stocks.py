import argparse
import pandas as pd
import yfinance as yf
from datetime import datetime
from dateutil.relativedelta import relativedelta
from typing import List

def import_ticker_list(input_file: str) -> List[str]:
    """Imports the list of stock tickers from the file provided.
    Filters out any empty lines. If the file could not be read or
    if no tickers could be found, returns None."""
    try:
        with open(input_file) as file:
            lines = file.read().splitlines()
            return list(set(filter(None, lines)))
    except:
        return None

def ticker_history(tickers: List[str], days=30) -> pd.DataFrame:
    """Downloads the price history of all stock tickers provided over
    the given duration. Returns the adjusted closing prices for all
    stocks over the duration in a Pandas DataFrame."""
    end = datetime.today().strftime("%Y-%m-%d")
    delta = relativedelta(days=days)
    start = (datetime.today()-delta).strftime("%Y-%m-%d")
    data = yf.download(tickers, start=start, end=end)
    df = pd.DataFrame(data["Adj Close"])
    # Remove rows with empty values
    df.dropna(inplace=True)
    return df

def main(input_file="dow.stocks", output_file="stock_history.csv", days=0):
    tickers = import_ticker_list(input_file)
    if not tickers:
        print("||| ERROR: ticker list could not be imported.")
        return
    df = ticker_history(tickers, days)
    df.to_csv(output_file)

# Parse arguments and run the main script
if __name__ == "__main__":
    parser = argparse.ArgumentParser(usage="%(prog)s [-h] [-i INPUT] [-o -OUTPUT] [-d DAYS]")
    parser.add_argument("-i", "--input", help="file path to stock tickers (default: dow.stocks)", type=str, default="dow.stocks", metavar="")
    parser.add_argument("-o", "--output", help="ouptut .csv file path (default: dow_history.csv)", type=str, default="dow_history.csv", metavar="")
    parser.add_argument("-d", "--days", help="days of stock history (default: 30)", type=int, default=30, metavar="")
    args = parser.parse_args()
    main(input_file=args.input, output_file=args.output, days=args.days)