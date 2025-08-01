# Binary Star Finder Instructions & Definitions

## 1. Parameter Definitions
---

- **RA Threshold (deg)**:  
  Maximum allowed difference in Right Ascension (in degrees) between stars to consider them as a pair.

- **Dec Threshold (deg)**:  
  Maximum allowed difference in Declination (in degrees) between stars for pairing.

- **RUWE Difference Threshold**:  
  Maximum allowed difference in RUWE (Renormalized Unit Weight Error) values. RUWE indicates data quality; similar RUWE values suggest similar data reliability. Larger RUWE values suggest a lower fit to the astrometric single-star solution, potentially indicating binary stars or other light sources. 

- **Distance Threshold (pc)**:  
  Maximum difference in distance (parsecs) between stars. Ensures physical proximity in space.

- **Minimum RUWE**:  
  Minimum RUWE value to filter stars with poor astrometric fits. Only stars with RUWE >= this value are considered.

## 2. Reasoning for Choosing Parameters
---

These thresholds help identify candidate binary stars by filtering for stars that are close together on the sky, at similar distances, and with similar data quality (RUWE). Tight thresholds reduce false positives but might miss wider binaries.

## 3. Gaia Catalog
---

Gaia is a space observatory that provides precise astrometric data (positions, parallaxes, proper motions) for over a billion stars. Your input CSV should contain at least:

- `source_id` (string): Unique star identifier.  
- `ra` (float): Right Ascension in degrees.  
- `dec` (float): Declination in degrees.  
- `parallax` (float): Parallax in milliarcseconds (mas), used to compute distance.  
- `ruwe` (float): Astrometric goodness-of-fit indicator.

**Example CSV header (minimum required columns):**
```
source_id,ra,dec,parallax,ruwe
1234567890123456789,150.1234,-35.5678,10.5,1.25
```

## 4. WDS Catalog
---

Washington Double Star Catalog (WDS) lists known binary stars with precise coordinates. The app compares Gaia stars to WDS entries to identify known binaries.

Your WDS catalog should be a text file with star coordinates at the end of each line, for example:

```
... 12 34 56.7 +12 34 56.7
```

Coordinates are in the format HH MM SS.S (RA) and ±DD MM SS.S (Dec).

## 5. File Formats
---

- **Gaia CSV**: Comma-separated values with headers, UTF-8 encoded, containing required columns.  
- **WDS TXT**: Plain text file with fixed-width lines ending with coordinate strings parseable by astropy's SkyCoord.

---

Feel free to adjust thresholds to balance sensitivity and precision. Contact Ansh Menghani at ansh.menghani@gmail.com for support.

## © 2025 Ansh Menghani. All rights reserved.
