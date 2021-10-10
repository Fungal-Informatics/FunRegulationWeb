import "../config/env-vars.ts";
import { Umzug } from "umzug";
import { getUmzugConfigured } from "../config/migration";

async function runCommands(w: string[], umz: Umzug) {
	// valid comands

	// run MIGRATION_NAME
	// run all
	// run many MIGRATION_NAME MIGRATION_NAME MIGRATION_NAME MIGRATION_NAME
	// run to MIGRATION_NAME
	// run from MIGRATION_NAME
	// run from MIGRATION_NAME to MIGRATION_NAME
	// revert last
	// revert all
	// revert MIGRATION_NAME
	// revert to MIGRATION_NAME
	// revert many MIGRATION_NAME MIGRATION_NAME MIGRATION_NAME MIGRATION_NAME
	const invalid = () => console.log("Invalid option");

	try {
		switch (w[0]) {
			case "run":
				switch (w[1]) {
					case "all":
						await umz.up();
						break;
					case "to":
						await umz.up({ to: w[2] });
						break;
					case "from":
						await umz.up({ from: w[2] });
						break;
					case "many":
						await umz.up(w.slice(2));
						break;
					default:
						invalid();
						break;
				}
				break;
			case "revert":
				switch (w[1]) {
					case "all":
						await umz.down({ to: 0 });
						break;
					case "last":
						await umz.down();
						break;
					case "to":
						await umz.down({ to: w[2] });
						break;
					case "many":
						await umz.down(w.slice(2));
						break;
					default:
						invalid();
						break;
				}
				break;
			default:
				invalid();
				break;
		}
	} catch (err) {
		console.log("ERROR: ", err);
	}

	/*
    run
        all
        many
            MIGRATION_NAME
                MIGRATION_NAME
                    MIGRATION_NAME
                        ...
        to
            MIGRATION_NAME
        from
            MIGRATION_NAME
        MIGRATION_NAME
    revert
        all
        last
        to
            MIGRATION_NAME
        many
            MIGRATION_NAME
                MIGRATION_NAME
                    MIGRATION_NAME
                        ...
        MIGRATION_NAME
     */
}

async function cliHelper(comands: string[]) {
	runCommands(comands.splice(2), await getUmzugConfigured());
}

cliHelper(process.argv);
