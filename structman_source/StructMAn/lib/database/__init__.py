from pkgutil import walk_packages


# import all modules' attributes into structman.lib.database namespace, could cause issues if modules have any overlapping attributes
for loader, module_name, _ in walk_packages(__path__):
    module = loader.find_module(module_name).load_module(module_name)
    module_dict = module.__dict__
    try:
        imports = module.__all__
    except AttributeError:
        imports = [attr for attr in module_dict if not attr.startswith('_')]
    globals().update({attr: module_dict[attr] for attr in imports})
