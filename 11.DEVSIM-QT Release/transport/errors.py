"""Common transport failure contract; SCF does not interpret backend diagnostics."""
class TransportDomainError(RuntimeError):
    status='MODEL_UNDERSPECIFIED'
    def __init__(self,message,diagnostics=None):
        self.diagnostics={} if diagnostics is None else diagnostics
        super().__init__(message)

class TransportNumericalError(RuntimeError):
    status='NUMERICAL_FAILURE'
