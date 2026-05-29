# data/

Catalog payload and provenance for the WCT-motivated GWTC diagnostic.

## Files

- `gwtc_allevents_raw.json` — verbatim payload from the official GWOSC Event
  API (`https://gwosc.org/eventapi/json/allevents/`). Regenerate with
  `python scripts/fetch_gwtc_catalog.py`.
- `manifest.csv` — provenance: source, URL, UTC access date, SHA-256, byte
  size, total event count, and per-catalog entry counts. Regenerated on every
  fetch.

## Source

Official **GWOSC / LIGO-Virgo-KAGRA** catalog metadata — not scraped secondary
summaries. See [../PROVENANCE.md](../PROVENANCE.md) for the full record.

## Notes

- The raw JSON is the single source of truth; `tables/gwtc_events_clean.csv`
  (built by `scripts/build_gwtc_table.py`) is derived from it.
- If GWOSC updates the upstream catalog, the SHA-256 in `manifest.csv` will
  change — re-run the fetch and record the new hash.
- `scripts/fetch_gwtc_catalog.py --offline-ok` reuses a cached payload when the
  network is unavailable.
