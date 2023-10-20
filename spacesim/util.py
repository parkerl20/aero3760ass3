class DictLogger(): 
    def __init__(self) -> None:
        self.log: dict = dict()
        
    def add_log(self, key: str, init_value: any) -> None:
        if key in self.log:
            return
        
        self.log[key] = init_value
        return
    
    def get_log(self, key: str) -> any:
        if key not in self.log:
            raise KeyError(f"Key {key} not in log")
        
        return self.log[key]
    
    def get_log_keys(self) -> list[str]:
        return list(self.log.keys())