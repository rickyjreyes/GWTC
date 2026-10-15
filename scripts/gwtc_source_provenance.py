#!/usr/bin/env python3
"""Explicit byte-only LF/CRLF audit for frozen CSVs historically hashed on Windows.

Never accept changed CSV field/event content. A match is accepted only if the
current file differs from the declared hash solely by uniform newline encoding;
the observed checkout SHA and exact transformed SHA are both reported.
"""
from __future__ import annotations
import hashlib
from pathlib import Path

def hashes(path):
    raw=Path(path).read_bytes()
    canonical=raw.replace(b"\r\n",b"\n")
    if b"\r" in canonical:
        raise ValueError("Nonstandard bare-CR bytes; manual source audit required")
    return {
      "checkout_sha256":hashlib.sha256(raw).hexdigest(),
      "normalized_lf_sha256":hashlib.sha256(canonical).hexdigest(),
      "normalized_crlf_sha256":hashlib.sha256(canonical.replace(b"\n",b"\r\n")).hexdigest(),
      "checkout_bytes":len(raw),
      "content_equivalence_scope":"Only LF-to-CRLF representation, not CSV parsing, field changes or line reordering"}

def validate_frozen_csv(path,expected):
    h=hashes(path)
    if expected==h["checkout_sha256"]:match="EXACT_ORIGINAL_BYTES"
    elif expected==h["normalized_crlf_sha256"] or expected==h["normalized_lf_sha256"]:
        match="DECLARED_HASH_MATCHES_ONLY_LF_CRLF_TRANSFORM"
    else:
        raise ValueError("Frozen CSV source mismatch beyond permitted newline transformation: "+str(path)+
            " expected="+str(expected)+" actual="+h["checkout_sha256"]+
            " LF="+h["normalized_lf_sha256"]+" CRLF="+h["normalized_crlf_sha256"])
    return {"file":str(path),"declared_original_sha256":expected,"verification_status":match,**h}
