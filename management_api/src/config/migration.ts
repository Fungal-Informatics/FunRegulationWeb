import { Umzug, migrationsList } from "umzug";
import { getMigrationList } from "../migrations/list";
import { DBMigrationStorage } from "./migration-db-storage";
import { pool } from "./pg-pool";

export async function getUmzugConfigured() {
	return new Umzug({
		migrations: migrationsList(await getMigrationList()),
		storage: new DBMigrationStorage(pool),
		logging: console.log,
	});
}

export async function runMigrations() {
	const umz = await getUmzugConfigured();
	await umz.up();
}
