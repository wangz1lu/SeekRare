#!/usr/bin/env python3
"""Standalone dbSNP common variant filter — converts dbSNP to match input VCF format."""

import argparse
from seekrare.preprocess.dbsnp_filter import run_dbsnp_filter

def main():
    p = argparse.ArgumentParser(description="dbSNP common 变异过滤")
    p.add_argument("--input", required=True, help="输入 VCF.gz")
    p.add_argument("--dbsnp", required=True, help="dbSNP common VCF.gz")
    p.add_argument("--output", required=True, help="输出 VCF.gz（剔除 common 后）")
    args = p.parse_args()
    r = run_dbsnp_filter(args.input, args.dbsnp, args.output)
    print(f"✅ {r['n_removed']}/{r['n_before']} removed (input={r['input_format']}, dbsnp={r['dbsnp_format']})")
    print(f"   Output: {args.output}")

if __name__ == "__main__":
    main()
