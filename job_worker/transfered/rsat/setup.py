"""
##################### RSAT WS ########################
"""
""" 
    Create the SOAP client to request RSAT services
""" 
# Define URL for RSAT services (alternative RSAT WSs to use when the UNAM WS is down)
#wsdlUrl =  'http://rsat.ulb.ac.be/rsat/web_services/RSATWS.wsdl'
#wsdlUrl =  'http://rsat01.biologie.ens.fr/rsat/web_services/RSATWS.wsdl'
#wsdlUrl = 'http://pedagogix-tagc.univ-mrs.fr/rsat/web_services/RSATWS.wsdl'
wsdlUrl = 'http://embnet.ccg.unam.mx/rsat/web_services/RSATWS.wsdl'
#wsdlUrl = 'http://rsat-tagc.univ-mrs.fr/rsat/web_services/RSATWS.wsdl'

# Create the client
try:
    client = Client(wsdlUrl)
except (Exception) as error:
    lib.log.info("Connection to RSAT Server failed!", error)
    raise
# Need the service interface to perform requests
rsat_service = client.service

"""
    Define client header (optional)
"""
userAgent = 'RSAT-Client/v%s (%s; Python %s; %s)' % (
    GRNInferenceVersion, 
    os.path.basename( __file__ ),
    platform.python_version(), 
    platform.system()
)
httpHeaders = {'User-agent': userAgent}
client.set_options(headers=httpHeaders)
client.set_options(timeout=300)

