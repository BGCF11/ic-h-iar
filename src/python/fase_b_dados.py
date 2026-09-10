"""
Phase B — Asynchronous extraction of the Sentinel-2 spatiotemporal grid
=====================================================================
This script implements the first step of the H-IAR Phase B workflow:
full extraction of the Mata de Santa Genebra ARIE spatiotemporal dataset
through Google Earth Engine.

The procedure defines a georeferenced study polygon, selects the harmonized
Sentinel-2 Level-2A Surface Reflectance collection, retains only pixels
classified as vegetation by the Scene Classification Layer (SCL == 4),
rescales the optical bands to percentage reflectance, and exports the full
table asynchronously to Google Drive.

References:
- Google Earth Engine Data Catalog: COPERNICUS/S2_SR_HARMONIZED.
- Google Earth Engine API: Image.sample, Image.pixelLonLat, Image.addBands,
  Image.updateMask and Export.table.toDrive.
"""

import ee

# Authenticate and initialize Google Earth Engine for the research project.
# This step enables submission of the asynchronous export task.
ee.Authenticate()
ee.Initialize(project='hiar-497604')

def exportar_malha_florestal_to_drive():
    """
    Submit an asynchronous Earth Engine export of the Sentinel-2 spatiotemporal
    dataset for the Mata de Santa Genebra ARIE.

    The function:
    1. defines the expanded study polygon, including the interior and edge;
    2. prepares each Sentinel-2 image with the SCL == 4 mask, optical bands,
       geographic coordinates, and a numerical timestamp band;
    3. samples pixels at a scale of 10 m;
    4. exports the resulting FeatureCollection to Google Drive as a CSV.

    Reference: GEE Export.table.toDrive; Sentinel-2 SR Harmonized.
    """
    # 1. Define the expanded Phase B polygon.
    # The vertices delimit the study area in the Mata de Santa Genebra ARIE,
    # including the forest interior and peripheral area for exploratory edge analysis.
    coordenadas_contorno = [
        [-47.119989, -22.809194], [-47.126428, -22.819750],
        [-47.123011, -22.820503], [-47.116847, -22.819675],
        [-47.115114, -22.821789], [-47.115861, -22.826542],
        [-47.108186, -22.833247], [-47.104822, -22.833867],
        [-47.104480, -22.823980], [-47.103386, -22.822903],
        [-47.1079777, -22.8165861]
    ]
    
    # Explicitly close the polygon ring:
    # make the last vertex coincide with the first before creating the geometry.
    if coordenadas_contorno[0] != coordenadas_contorno[-1]:
        coordenadas_contorno.append(coordenadas_contorno[0])

    # Polygon geometry corresponding to the georeferenced boundary.
    limites_mata = ee.Geometry.Polygon([coordenadas_contorno])

    def processar_imagem(imagem):
        """
        Prepare a Sentinel-2 image for tabular sampling.

        Retain vegetation pixels using SCL, select the four optical bands used
        as H-IAR components, add numerical longitude/latitude bands, and insert
        system:time_start as the constant band 'time'.

        Representing time as a band is necessary because img.sample() exports
        properties by band; each sampled row therefore retains the acquisition
        time of its source image.
        """
        # SCL mask: class 4 = Vegetation in the Sentinel-2 L2A catalogue.
        # Retain only pixels classified as vegetation.
        mascara_vegetacao = imagem.select('SCL').eq(4) 

        # Spectral components of the H-IAR vector:
        # B2 = blue, B3 = green, B4 = red, B8 = near infrared.
        #
        # Sentinel-2 SR bands in GEE are stored as DN with a scale factor of 0.0001
        # for physical reflectance. Dividing by 100 converts DN to percentage
        # reflectance: DN/100 = 100 * physical reflectance.
        #
        # Percentage scaling keeps the numerical magnitude suitable for the
        # covariance matrices subsequently used by the Kalman filter.
        bandas = imagem.select(['B2', 'B3', 'B4', 'B8']).divide(100)

        # Pixel coordinates, incorporated as numerical bands.
        coordenadas_grid = ee.Image.pixelLonLat()
        
        # Insert the raw system:time_start timestamp as a constant band named 'time'.
        # Since img.sample() exports one property per image band, this approach
        # ensures that each sampled pixel also carries the scene acquisition time.
        banda_tempo = ee.Image.constant(imagem.getNumber('system:time_start')).rename('time')
        
        # Combine reflectance, coordinates and time in a single multiband image.
        imagem_preparada = bandas.addBands(coordenadas_grid).addBands(banda_tempo)

        # Apply the vegetation mask to the multiband image.
        # Pixels outside SCL class 4 are masked and subsequently removed
        # by dropNulls=True in sample().
        return imagem_preparada.updateMask(mascara_vegetacao)

    # 2. Request the Sentinel-2 SR Harmonized collection.
    # The harmonized collection corrects the radiometric offset associated with
    # PROCESSING_BASELINE >= 04.00, aligning older and newer scenes.
    #
    # filterBounds restricts the collection to the study polygon;
    # filterDate sets [2020-01-01, 2023-12-31), with an exclusive end date;
    # map applies the reflectance and table preparation to each scene.
    colecao = (ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
               .filterBounds(limites_mata)
               .filterDate('2020-01-01', '2023-12-31')
               .map(processar_imagem))

    # 3. Spatial sampling at 10 m.
    # For each Sentinel-2 scene, img.sample() converts valid pixels within the
    # polygon into a tabular FeatureCollection with one property per band:
    # longitude, latitude, B2, B3, B4, B8 and time.
    #
    # map applies this sampling to all scenes; flatten concatenates the
    # individual FeatureCollections into one spatiotemporal table.
    amostras_tabelares = colecao.map(lambda img: img.sample(
        region=limites_mata,
        scale=10,
        projection='EPSG:4326',
        geometries=False,
        dropNulls=True
    )).flatten()

    # 4. Create the asynchronous export task.
    # Processing is delegated to the GEE batch task system,
    # avoiding the practical limits of synchronous local requests.
    print("A preparar a tarefa de exportação. Este processo pode demorar alguns segundos...")
    
    tarefa = ee.batch.Export.table.toDrive(
        collection=amostras_tabelares,
        description='Extracao_MataSantaGenebra_Integral_Corrigida',
        folder='IC_Sensoriamento_H_IAR',
        fileFormat='CSV',
        # selectors fixes the order and subset of exported columns.
        # 'time' is included because processar_imagem() added it as a constant band.
        selectors=['longitude', 'latitude', 'B2', 'B3', 'B4', 'B8', 'time'] 
    )
    
    # Execution now proceeds asynchronously in Earth Engine.
    tarefa.start()
    print("Tarefa submetida com sucesso aos servidores da Google Cloud!")

# Direct execution of the Phase B extraction stage.
exportar_malha_florestal_to_drive()
