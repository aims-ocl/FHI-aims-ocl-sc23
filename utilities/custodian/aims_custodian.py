from custodian import Custodian
from custodian.fhi_aims.jobs import AimsJob
from custodian.fhi_aims.handlers import ConvergenceEnhancer, AimsRelaxHandler, FrozenJobErrorHandler
from custodian.fhi_aims.validators import AimsConvergedValidator

job = [AimsJob(aims_cmd=['mpirun', 'aims'])]
handlers = [AimsRelaxHandler(), FrozenJobErrorHandler()]  # , ConvergenceEnhancer()]
valids = [AimsConvergedValidator()]

cus = Custodian(handlers=handlers, validators=valids, jobs=job)

cus.run()

