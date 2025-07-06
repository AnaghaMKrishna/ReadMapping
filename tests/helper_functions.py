import pytest

def check_class_method(cls, meth):
    """Check class for a method and print a clear message if absent

    Args:
        meth (str): The method of which to assess the existence callable-ity
    """
    __tracebackhide__ = True # don't show this function in the traceback so just the informative test is shown
    # Figure out if message should say a or an for this method
    a_an = "an" if meth[0] in {"a", "e", "i", "o" "u"} else "a"
    if not hasattr(cls, meth):
        pytest.fail(f"Your {cls.__name__} class does not have {a_an} {meth} method.")
    else:
        if not callable(eval(f"cls.{meth}")):
            pytest.fail(f"Your {cls.__name__} class does have {a_an} {meth} attribute, but it is not callable")

def check_class_attribute(cls, attr):
    """Check class for an attribute and print a clear message if absent

    Args:
        attr (str): The attribute of which to assess the existence
    """
    __tracebackhide__ = True
    a_an = "an" if attr[0] in {"a", "e", "i", "o" "u"} else "a"
    if not hasattr(cls, attr):
        pytest.fail(f"Your {cls.__name__} class does not have {a_an} {attr} attribute.")

def check_attribute_value(instance, attr, expected, tested_data):
    """Compare an attribute of a class instance to an expected value

    Args:
        instance (Object): The instance you want to assess
        attr (str): The attribute whose value should be compared
        expected (any): The expected value
        tested_data (str): the nature of the tested data
    """
    __tracebackhide__ = True
    a_an = "an" if tested_data[0] in {"a", "e", "i", "o", "u"} else "a"
    value = getattr(instance, attr)
    try:
        assert value == expected
    except:
        pytest.fail(f"Your {instance.__class__.__name__}.{attr} contains {value} for {a_an} {tested_data}, when it should have contained {expected}.")


def check_method_output(instance, method, expected, tested_data, args=(), kwargs=None):
    """Compare an attribute of a class instance to an expected value

    Args:
        instance (Object): The instance you want to assess
        method (str): The method whose return value should be compared
        args (tuple[any]): The arguments to provide to the method when called
        kwargs (dict[str, any]): The keyword arguments to provide to the method when called
        expected (any): The expected value
        tested_data (str): the nature of the tested data
    """
    __tracebackhide__ = True
    if kwargs == None:
        kwargs = {}
    a_an = "an" if tested_data[0] in {"a", "e", "i", "o", "u"} else "a"
    try:
        value = eval(f"instance.{method}(*args, **kwargs)")
    except Exception as e:
        pytest.fail(f"Your {instance.__class__.__name__}.{method} for {a_an} {tested_data}, when it was run with input {args}. The error was {e}")
    try:
        assert value == expected
    except:
        pytest.fail(f"Your {instance.__class__.__name__}.{method} returned {repr(value)} for {a_an} {tested_data}, when it should have returned {repr(expected)}.")
