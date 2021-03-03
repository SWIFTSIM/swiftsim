#!/usr/bin/env python3
import sys
import re
import typing
from task import tasks as ref_tasks


def get_string_between(s: str, beg: str, end: str):
    """
    Match the smallest string contained within
    the substring beg and end. Returns None if nothing
    is found
    """
    regex = beg + ".*?" + end
    m = re.search(regex, s)

    if m is None:
        return None

    m = m.group(0)
    return m[len(beg):-len(end)]


class Task:
    """
    Class representing a single task along with its links.
    """
    def __init__(self, s: str):
        """
        Initialize the task from a line in the input file.
        """

        self.name = s[:s.find("[")]
        self.iact_type = self.read_attribute("type", s)
        self.active = self.read_attribute("active", s)
        self.level = self.read_attribute("level", s)
        self.updates = self.read_attribute("updates", s)
        self.policy = self.read_attribute("policy", s)
        self.cell_contains = self.read_attribute("cell_contains", s)
        self.act_on = self.read_attribute("act_on", s)
        self.link_to = []
        self.linked_from = []

    def read_attribute(self, attribute_name: str, s: str):
        """
        Read a single attribute in the definition line.
        """
        attribute_name += "="
        attr = get_string_between(s, attribute_name, ",")
        if attr is not None:
            return attr

        attr = get_string_between(s, attribute_name, "]")
        return attr

    def add_linked_from(self, name: str, task_type: str):
        """
        Add a link for this task.
        """
        self.linked_from.append({"name": name,
                                 "task_type": task_type})

    def add_link_to(self, name: str, task_type: str):
        """
        Add a link for this task.
        """
        self.link_to.append({"name": name,
                             "task_type": task_type})

    def print_task(self):
        """
        Print the task.
        """
        print("{}: iact={}, active={}, level={}, link={}".format(
            self.name, self.iact_type, self.active,
            self.level, self.link_to))

    def write_if_condition(self, f: typing.TextIO, level: bool = True,
                           policy: bool = True, task_type: bool = False,
                           task_subtype: bool = False):
        """
        Write the if condition for this task.
        The conditions can be disabled if needed.

        Returns
        -------

        written: bool
            Did this function write something (e.g. empty if condition)?
        """
        if_condition = ""
        if self.level and level:
            if_condition += "(c->%s == c) && " % self.level
        if self.policy and policy:
            if_condition += "(e->policy & engine_policy_%s) && " % self.policy
        if task_type:
            task_type_name = task_type
            if isinstance(task_type, bool):
                task_type_name = ref_tasks[self.name]["type"]
            if_condition += "(t->type == task_type_%s) && " % \
                task_type_name
        if task_subtype:
            task_subtype_name = task_subtype
            if isinstance(task_subtype, bool):
                task_subtype_name = ref_tasks[self.name]["subtype"]

            if_condition += "(t->subtype == task_subtype_%s) && " % \
                task_subtype_name

        if len(if_condition) == 0:
            return False

        # Remove ' &&'
        if_condition = if_condition[:-4]

        if_condition = "if (%s) {\n" % if_condition
        f.write(if_condition)

        return True

    def write_maketask_definitions(self, f: typing.TextIO):
        """
        Write the code corresponding to the definition of the task
        (e.g. scheduler_addtask).
        It creates the tasks acting on a single cell at a time and
        does not care about dependencies.
        This will replace the work done in engine_make_hierarchical_tasks.
        """
        is_implicit = int(self.iact_type == "implicit")
        if self.iact_type != "single" and not is_implicit:
            return

        # Compute the scheduler_addtask
        task_type = ref_tasks[self.name]["type"]
        creation = "c->{name} = scheduler_addtask(s, task_type_{task_type},"
        creation += " task_subtype_none, 0, {implicit}, c, NULL)"
        creation = creation.format(
            name=ref_tasks[self.name]["variable"], task_type=task_type,
            implicit=is_implicit)

        # Write the code inside a file
        if_done = self.write_if_condition(f)
        code = "\t %s;\n" % creation
        if if_done:
            code += "};\n"
        f.write(code)

    def write_maketask_recursions(
            self, f: typing.TextIO, tasks: dict):
        if self.iact_type != "splitted":
            return

        func = "void engine_add_{name}(struct engine *e, struct cell *c,\n"
        func += "\tstruct task *task_in, struct task *task_out) {{\n\n"
        func += "\tif (c->hydro.count_total == 0) return;\n"
        func += "\tif (!c->split || c->hydro.count_total < engine_max_parts_per_{name}) {{\n\n"
        func += "\t\tstruct scheduler *s = &e->sched;\n"
        func += "\t\tc->{variable} = scheduler_addtask(s, task_type_{task_type}, task_subtype_none,\n"
        func += "\t\t\t0, 0, c, NULL);\n"
        func += "\t\tscheduler_addunlock(s, task_in, c->{variable});\n"
        func += "\t\tscheduler_addunlock(s, c->{variable}, task_out);\n"
        func += "\t}} else {{\n"
        func += "\t\tfor(int k = 0; k < 8; k++)\n"
        func += "\t\t\tif(c->progeny[k] != NULL)\n"
        func += "\t\t\t\tengine_add_{name}(e, c->progeny[k], task_in, task_out);\n"
        func += "\t}}\n}}\n\n"
        func = func.format(
            name=self.name, variable=ref_tasks[self.name]["variable"],
            task_type=ref_tasks[self.name]["type"])
        f.write(func)

    def write_maketask_single_single(
            self, f: typing.TextIO, link: str, tasks: dict):
        # Generate part of the variable name
        self_level = ""
        if self.level:
            self_level = "%s->" % self.level

        # Grab data for the link
        link = tasks[link]
        link_level = ""
        if link.level:
            link_level = "%s->" % link.level

        # Write the if condition
        write_policy = self.policy != link.policy
        if_done2 = link.write_if_condition(
            f, level=False, policy=write_policy)

        # Generate the scheduler_addunlock
        unlock = "\t\tscheduler_addunlock(s, ci->{level1}{name1},"
        unlock += " ci->{level2}{name2});\n"
        unlock = unlock.format(
            level1=self_level, name1=ref_tasks[self.name]["variable"],
            level2=link_level, name2=ref_tasks[link.name]["variable"]
        )

        # Close parenthesis
        if if_done2:
            unlock += "\t};\n"
        f.write(unlock)

    def write_splitted(
            self, f: typing.TextIO, tasks: dict):

        if len(self.link_to) != 1 or len(self.linked_from) != 1:
            raise Exception("Cannot recurse with more than 1 link")

        t_bef = tasks[self.linked_from[0]["name"]]
        t_aft = tasks[self.link_to[0]["name"]]

        level = ""
        if t_bef.level:
            level = t_bef.level + "->"
        before = "ci->{level}{name}".format(
            level=level, name=ref_tasks[t_bef.name]["variable"])

        level = ""
        if t_aft.level:
            level = t_aft.level + "->"
        after = "ci->{level}{name}".format(
            level=level, name=ref_tasks[t_aft.name]["variable"])


        # Write the task
        if_written = self.write_if_condition(
            f, level=False, policy=True, task_type=ref_tasks[t_bef.name]["type"])

        add = "\tengine_add_{name}(e, ci, {before}, {after});\n"
        f.write(add.format(
            name=self.name, before=before, after=after))

        if if_written:
            f.write("}\n")

    def write_maketask_extra_loop_single(
            self, f: typing.TextIO, tasks: dict):
        """
        Write the code corresponding to link_to for a single cell task.
        """
        if len(self.link_to) == 0:
            return

        iact1 = self.iact_type
        if iact1 == "iact":
            raise Exception("Should not happen")

        # Write the initial if condition
        if_done = self.write_if_condition(
            f, level=False, task_type=True)

        # loop over all the links
        for link in self.link_to:
            link_name = link["name"]

            iact2 = tasks[link_name].iact_type

            # Only link tasks that are not iact.
            if iact2 == "iact" or iact2 == "splitted":
                continue

            self.write_maketask_single_single(
                f, link_name, tasks)

        # Close parenthesis
        if if_done:
            f.write("};\n")

    def write_iact_link_pair(self, f, tasks, task_type):
        if_written = self.write_if_condition(
            f, level=False, policy=True, task_type=task_type, task_subtype="density")

        # Create the task
        if self.name != "density":
            # Create the task
            create = "\tstruct task *new = scheduler_addtask(s,"
            create += " task_type_{task_type}, task_subtype_{subtype}, "
            create += "flags, 0, ci, cj);\n"
            f.write(create.format(task_type=task_type,
                                  subtype=ref_tasks[self.name]["subtype"]))

            # Add it to the link
            addlink = "\tengine_addlink(e, &{cell}->{variable}, new);\n"
            f.write(addlink.format(
                cell="ci", variable=ref_tasks[self.name]["variable"]))
            f.write(addlink.format(
                cell="cj", variable=ref_tasks[self.name]["variable"]))

        else:
            f.write("\tstruct task *new = t;\n")

        # Make link towards current task
        for link_from in self.linked_from:
            link_name = link_from["name"]

            if link_from["task_type"] is not None and link_from["task_type"] != task_type:
                continue

            t2 = tasks[link_name]
            iact2 = t2.iact_type
            if iact2 == "iact":
                raise Exception("Should not happen")

            link_level = ""
            if t2.level:
                link_level = "%s->" % t2.level

            unlock = "\tscheduler_addunlock(s, {cell}->{level}{name}, new);\n"
            f.write(unlock.format(
                cell="ci", level=link_level,
                name=ref_tasks[link_name]["variable"]))

            if t2.level:
                f.write("\tif(ci->%s != cj->%s) {\n" % (t2.level, t2.level))

            f.write("\t" + unlock.format(
                cell="cj", level=link_level,
                name=ref_tasks[link_name]["variable"]))

            if t2.level:
                f.write("}\n")

        # Make link from current task
        for link_to in self.link_to:
            link_name = link_to["name"]

            if link_to["task_type"] is not None and link_to["task_type"] != task_type:
                continue

            t2 = tasks[link_name]
            iact2 = t2.iact_type
            if iact2 == "iact":
                raise Exception("Should not happen")

            link_level = ""
            if t2.level:
                link_level = "%s->" % t2.level

            unlock = "\tscheduler_addunlock(s, new, {cell}->{level}{name});\n"
            f.write(unlock.format(
                cell="ci", level=link_level,
                name=ref_tasks[link_name]["variable"]))

            if t2.level:
                f.write("\tif(ci->%s != cj->%s) {\n" % (t2.level, t2.level))

            f.write("\t" + unlock.format(
                cell="cj", level=link_level,
                name=ref_tasks[link_name]["variable"]))

            if t2.level:
                f.write("}\n")

        if if_written:
            f.write("}\n")

    def write_iact_link_self(self, f, tasks, task_type):
        if_written = self.write_if_condition(
            f, level=False, policy=True, task_type=task_type, task_subtype="density")

        # Create the task
        if self.name != "density":
            # Create the task
            create = "\tstruct task *new = scheduler_addtask(s,"
            create += " task_type_{task_type}, task_subtype_{subtype}, "
            create += "flags, 0, ci, NULL);\n"
            f.write(create.format(task_type=task_type,
                                  subtype=ref_tasks[self.name]["subtype"]))

            # Add it to the link
            addlink = "\tengine_addlink(e, &ci->{variable}, new);\n"
            f.write(addlink.format(variable=ref_tasks[self.name]["variable"]))

        else:
            f.write("\tstruct task *new = t;\n")

        # Make link towards current task
        for link_from in self.linked_from:
            link_name = link_from["name"]

            if link_from["task_type"] is not None and link_from["task_type"] != task_type:
                continue

            t2 = tasks[link_name]
            iact2 = t2.iact_type
            if iact2 == "iact":
                raise Exception("Should not happen")

            link_level = ""
            if t2.level:
                link_level = "%s->" % t2.level

            unlock = "\tscheduler_addunlock(s, ci->{level}{name}, new);\n"
            f.write(unlock.format(level=link_level,
                                  name=ref_tasks[link_name]["variable"]))

        # Make link from current task
        for link_to in self.link_to:
            link_name = link_to["name"]

            if link_to["task_type"] is not None and link_to["task_type"] != task_type:
                continue

            t2 = tasks[link_name]
            iact2 = t2.iact_type
            if iact2 == "iact":
                raise Exception("Should not happen")

            link_level = ""
            if t2.level:
                link_level = "%s->" % t2.level

            unlock = "\tscheduler_addunlock(s, new, ci->{level}{name});\n"
            f.write(unlock.format(level=link_level,
                                  name=ref_tasks[link_name]["variable"]))
        if if_written:
            f.write("}\n")

    def write_maketask_extra_loop_iact(
            self, f: typing.TextIO, tasks: dict):
        if len(self.link_to) == 0:
            raise Exception("Should not happen")

        # Self case
        f.write("// self case\n")
        self.write_iact_link_self(f, tasks, task_type="self")

        # Subself case
        f.write("// sub self case\n")
        self.write_iact_link_self(f, tasks, task_type="sub_self")

        # Pair case
        f.write("// pair case\n")
        self.write_iact_link_pair(f, tasks, task_type="pair")

        # Subpair case
        f.write("// sub pair case\n")
        self.write_iact_link_pair(f, tasks, task_type="sub_pair")

    def write_maketask_extra_loop(
            self, f: typing.TextIO, tasks: dict):
        """
        Write the code corresponding to the link from this task
        (e.g. scheduler_addunlock)
        and defines it if it is a task acting on two cells.
        In order to have enough information for the links,
        a dictionary containing all the tasks is required.
        This will replace the work done in
        engine_make_extra_hydroloop_tasks.
        """

        # Write header
        f.write("// %s\n" % self.name)

        if self.iact_type == "iact":
            self.write_maketask_extra_loop_iact(f, tasks)
        elif self.iact_type != "splitted":
            self.write_maketask_extra_loop_single(f, tasks)
        else:
            self.write_splitted(f, tasks)


class Reader:
    """
    Class dealing with the task system.
    """
    def __init__(self, filename: str):
        """
        Initialize the reader from a .task file
        """
        self.filename = filename
        self.tasks = {}
        self.name = None
        self.read()

    def read(self):
        """
        Read the .task file
        """
        with open(self.filename, "r") as f:
            for line in f:
                line = line.rstrip()
                self.read_line(line)

    def read_line(self, line: str):
        """
        Read a single line in the .task file
        """
        if "->" in line:
            self.read_link(line)
        elif "[" in line and "]" in line:
            self.read_definition(line)
        elif "label" in line:
            s = line.find('"')
            e = line.rfind('"')
            self.name = line[s+1:e]

    def read_link(self, line: str):
        """
        Read a line containing a link
        """
        names = line.split("->")
        if names[0] not in self.tasks:
            raise Exception(
                "Trying to link {} without any definition".format(
                    names[0]))

        task_type = None

        # extract the task name
        if "pair_" in names[0] or "self_" in names[0]:
            split = names[0].split("_")
            names[0] = split[-1]
            task_type = "_".join(split[:-1])
        if "pair_" in names[1] or "self_" in names[1]:
            if task_type is not None:
                raise Exception("Not implemented")
            split = names[1].split("_")
            names[1] = split[-1]
            task_type = "_".join(split[:-1])

        self.tasks[names[0]].add_link_to(names[1], task_type=task_type)
        self.tasks[names[1]].add_linked_from(names[0], task_type=task_type)

    def read_definition(self, line: str):
        """
        Read a line containing a definition.
        """
        t = Task(line)
        self.tasks[t.name] = t

    def print_reader(self):
        """
        Print the task system
        """
        print("Simulation type:", self.name)
        for t in self.tasks:
            self.tasks[t].print_task()

    def write_maketask(self, filename: str):
        """
        Write the code corresponding to the task system into a file.
        """
        with open(filename, "w") as f:
            f.write("// Hierarchical taks\n")
            for name in self.tasks:
                self.tasks[name].write_maketask_definitions(f)

            f.write("\n")
            f.write("// Recursions\n")
            for name in self.tasks:
                self.tasks[name].write_maketask_recursions(f, self.tasks)

            f.write("\n")
            f.write("// Dependencies\n")
            for name in self.tasks:
                self.tasks[name].write_maketask_extra_loop(f, self.tasks)


if __name__ == "__main__":
    reader = Reader(sys.argv[-1])

    reader.write_maketask("test.c")