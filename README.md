<!-- Made from the template weichseltree/atlas-layer. tools/atlas.py is an unchanged copy of lib/atlas.py in weichseltree/marbles. -->
# atlas-earth

Three layers of the real Earth for the atlas on weichseltree.com: **land**, **residents** and **height**, on
one latitude-longitude grid of 2048 × 4096 cells (about 9.8 km at the equator), as a tile pyramid of levels
0 to 3. The page reads `atlas.json` from this repository's site and draws the numbers; nothing of this
repository runs on the page.

| layer | each cell holds | source | licence |
|---|---|---|---|
| land | the share of its area that is not ocean (lakes and ice shelves count as land) | OpenStreetMap land polygons of 2026-01-03, via osmdata.openstreetmap.de | ODbL 1.0 |
| people | residents per km², 2025 | GHS-POP R2023A, epoch 2025, 30″ (European Commission, JRC) | CC BY 4.0 |
| height | mean height above sea level in metres, the sea floor negative | ETOPO 2022, 60″ surface (NOAA NCEI) | public domain |

**Attribution, as the sources ask:**

- Land: © OpenStreetMap contributors. The land layer is a database derived from OpenStreetMap and is
  available under the Open Database License (ODbL) 1.0: https://www.openstreetmap.org/copyright
- Residents: Schiavina, M., Freire, S., Carioli, A., MacManus, K. (2023): GHS-POP R2023A - GHS population
  grid multitemporal (1975-2030). European Commission, Joint Research Centre (JRC).
  doi:10.2905/2FF68A52-5B5B-4A22-8F40-C41DA8332CFE
- Heights: NOAA National Centers for Environmental Information (2022): ETOPO 2022 15 Arc-Second Global Relief
  Model (here its 60″ surface grid). doi:10.25921/fd45-gt74

The code of this repository (tools/, tests/) is under the MIT licence in `LICENSE`; the data is under the
licences above.

## How the arrays were made

`tools/prepare.py` made `results/earth/*.npy` once, on the machine that holds the sources, and wrote
`results/earth/summary.json` with the sources' sha256 and the totals it checked (the world's residents
kept to the person, the land area, the highest and lowest cell). It needs rasterio, netCDF4, pyogrio and
shapely and is never run by the Action.

- Land: the polygons drawn onto a Mercator grid of 16384 × 16384 samples; each cell takes the share of
  samples whose centres fall into it. North of 85.05° there is only sea, south of it only Antarctica.
- Residents: each 30″ cell's count goes to the cell holding its centre, then divided by the cell's area.
- Height: the mean of the 60″ cells whose centres fall into the cell.

How the stored numbers are read: a stored value `v` is `v * scale + offset` in the layer's unit (the
manifest gives both). Residents are stored in 65 536 steps from 0 to the densest cell (85 122 per km²),
a step of 1.3 per km², so a cell with fewer than about 0.65 residents per km² is stored as 0.

## Publishing

Push to `main`: the Action builds `dist/` from `results/` and `layers.toml`, checks it, runs the tests and
publishes it on this repository's Pages site. See the template's README for the details of the format.
