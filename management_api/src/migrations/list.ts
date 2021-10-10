import { MigrationDefinitionWithName } from "umzug/lib/migrationsList";

export async function getMigrationList(): Promise<
	MigrationDefinitionWithName[]
> {
	return [
		(await import("./00-users")).default,
		(await import("./01-boards")).default,
		(await import("./02-board_sharings")).default,
		(await import("./03-tasks")).default,
		(await import("./04-notes")).default,
		(await import("./05-bookmarks")).default,
		(await import("./06-sessions")).default,
	];
}
