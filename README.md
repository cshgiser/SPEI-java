# SPEI Calculator (Java)

A Java code implementation for calculating the Standardized Precipitation Evapotranspiration Index (SPEI) from climate data.

## Usage

To calculate SPEI:
1. Prepare your input data as 2D arrays:
   - `pre`: Precipitation data structured as `[years][months]`
   - `temp`: Temperature data structured as `[years][months]`
2. Create an SPEII instance with:
   - Your climate data arrays
   - Latitude of the location
   - Time scale (e.g., 1 for 1-month SPEI)
   - Starting year of your data
3. Retrieve results using the `getResult()` method

Example parameters:
- Latitude (`lat`): Location's latitude in degrees
- Scale: 1 (for 1-month SPEI)
- Start year: 1981

## About SPEI

The **Standardized Precipitation Evapotranspiration Index (SPEI)** is a multi-scalar drought index that combines precipitation and temperature data to assess wet/dry conditions. Key features:

- Calculated across various time scales (1, 3, 12, 24 months, etc.)
- Higher values indicate wetter conditions, lower values indicate drought
- Improves upon SPI by incorporating temperature effects
- Particularly suitable for climate change studies and regions with significant warming trends

### Reference
Vicente-Serrano, S. M., Beguería, S., & López-Moreno, J. I. (2010). A Multiscalar Drought Index Sensitive to Global Warming: The Standardized Precipitation Evapotranspiration Index. *Journal of Climate*, 23(7), 1696–1718. [DOI](https://doi.org/10.1175/2009JCLI2909.1)
